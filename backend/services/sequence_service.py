"""
Sequence Service - Fetch RefSeq track from UCSC,
build CDS using exon coordinates AND phases.
"""

import httpx
from typing import Optional, Dict, List
from loguru import logger
from tenacity import retry, stop_after_attempt, wait_exponential


class SequenceService:
    def __init__(self):
        self.base_url = "https://api.genome.ucsc.edu"
        self.cache: Dict[str, Dict] = {}

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=5))
    async def _fetch_track(self, chrom: str, start: int, end: int,
                           genome: str = "hg38") -> Optional[list]:
        """Fetch RefSeq track records from UCSC."""
        url = (f"{self.base_url}/getData/track"
               f"?genome={genome};track=refSeqComposite;chrom={chrom};start={start};end={end}")
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(url)
            r.raise_for_status()
            data = r.json()
            return data.get("refSeqComposite", [])

    def _pick_mane(self, records: list, mane_id: str) -> Optional[dict]:
        """Pick the MANE Select transcript from UCSC records."""
        for rec in records:
            if rec.get("name") == mane_id:
                return rec
        return None

    async def get_cds(self, gene: str, mane_id: str, chrom: str,
                      tx_start: int, tx_end: int,
                      genome: str = "hg38") -> Optional[Dict]:
        """
        Fetch UCSC RefSeq record for MANE transcript,
        assemble CDS using exonFrames.
        """
        if gene in self.cache:
            return self.cache[gene]

        # 1. Fetch the RefSeq track
        records = await self._fetch_track(chrom, tx_start - 100, tx_end + 100, genome)
        if not records:
            logger.error(f"❌ No RefSeq records for {chrom}:{tx_start}-{tx_end}")
            return None

        # 2. Pick MANE Select
        mane = self._pick_mane(records, mane_id)
        if not mane:
            logger.warning(f"⚠️ MANE {mane_id} not found, using first record")
            mane = records[0]

        logger.info(f"✅ Using UCSC RefSeq: {mane['name']}")

        # 3. Parse
        cds_start = mane["cdsStart"]        # 0-based
        cds_end = mane["cdsEnd"]            # 0-based, half-open
        strand = -1 if mane["strand"] == "-" else 1
        exon_starts = [int(x) for x in mane["exonStarts"].rstrip(",").split(",")]
        exon_ends = [int(x) for x in mane["exonEnds"].rstrip(",").split(",")]
        exon_frames = [int(x) for x in mane["exonFrames"].rstrip(",").split(",")]

        # 4. Build CDS by iterating exons in transcript order
        n = len(exon_starts)
        order = range(n - 1, -1, -1) if strand == -1 else range(n)

        fragments: List[str] = []
        for i in order:
            e_start = exon_starts[i]
            e_end = exon_ends[i]
            frame = exon_frames[i]

            # Skip non-coding exons
            if frame == -1:
                continue

            # Clip exon to CDS bounds
            frag_start = max(e_start, cds_start)
            frag_end = min(e_end, cds_end)
            if frag_start >= frag_end:
                continue

            # Fetch sequence from UCSC (0-based, half-open)
            seq = await self._fetch_dna(chrom, frag_start, frag_end, genome)
            if not seq:
                logger.error(f"❌ DNA fetch failed at {chrom}:{frag_start}-{frag_end}")
                return None

            # Apply frame offset (skip first `frame` bases if phase > 0)
            # For minus strand, the phase applies to the *reverse-complemented* fragment
            if strand == -1:
                seq = self._revcomp(seq)
                if frame > 0:
                    seq = seq[frame:]
            else:
                if frame > 0:
                    seq = seq[frame:]

            fragments.append(seq)

        cds = "".join(fragments)
        protein = self._translate(cds)

        logger.info(f"✅ {gene}: CDS={len(cds)} bp (expected 5592), protein={len(protein)} aa (expected 1863)")

        result = {
            "cds": cds,
            "protein": protein,
            "transcript_id": mane["name"],
            "gene": gene,
            "strand": strand,
            "source": "ucsc_refseq_mane",
        }
        self.cache[gene] = result
        return result

    async def _fetch_dna(self, chrom: str, start: int, end: int, genome: str) -> Optional[str]:
        url = (f"{self.base_url}/getData/sequence"
               f"?genome={genome};chrom={chrom};start={start};end={end}")
        try:
            async with httpx.AsyncClient(timeout=20.0) as client:
                r = await client.get(url)
                r.raise_for_status()
                seq = r.json().get("dna")
                return seq.upper() if seq else None
        except Exception as e:
            logger.error(f"❌ UCSC DNA fetch failed: {e}")
            return None

    @staticmethod
    def _revcomp(seq: str) -> str:
        return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def _translate(self, cds: str) -> str:
        codon_table = {
            'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T',
            'ACG':'T','ACT':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K',
            'AGC':'S','AGT':'S','AGA':'R','AGG':'R','CTA':'L','CTC':'L',
            'CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
            'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R',
            'CGG':'R','CGT':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V',
            'GCA':'A','GCC':'A','GCG':'A','GCT':'A','GAC':'D','GAT':'D',
            'GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
            'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F',
            'TTA':'L','TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
            'TGC':'C','TGT':'C','TGA':'*','TGG':'W',
        }
        protein = []
        for i in range(0, len(cds) - 2, 3):
            codon = cds[i:i+3]
            aa = codon_table.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        return ''.join(protein)


sequence_service = SequenceService()