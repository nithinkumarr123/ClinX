"""
Ensembl Service - Exon metadata from Ensembl, DNA from UCSC, Protein from Ensembl.
"""

import httpx
from typing import Optional, Dict
from loguru import logger
from config import config
from services.ucsc_service import ucsc_service


FALLBACK_CDS = {
    # ... keep your hardcoded BRCA1 / TP53 CDS + protein ...
}


class EnsemblService:
    def __init__(self):
        self.base_url = config.ENSEMBL_API_URL
        self.cds_cache = {}
        self.known_gene_ids = {
            "BRCA1": "ENSG00000012048",
            "TP53": "ENSG00000141510",
        }

    async def _fetch_transcript_metadata(self, gene: str) -> Optional[Dict]:
        """Fetch transcript + exons from Ensembl."""
        try:
            async with httpx.AsyncClient(timeout=30.0) as client:
                headers = {"Content-Type": "application/json"}

                r = await client.get(
                    f"{self.base_url}/lookup/symbol/homo_sapiens/{gene}",
                    headers=headers,
                )
                if r.status_code != 200:
                    return None
                gene_id = r.json().get("id")
                if not gene_id:
                    return None

                r = await client.get(
                    f"{self.base_url}/lookup/id/{gene_id}?expand=1",
                    headers=headers,
                )
                r.raise_for_status()
                data = r.json()
                transcripts = data.get("Transcript", [])
                if not transcripts:
                    return None

                canonical = next((t for t in transcripts if t.get("is_canonical")), None)
                if not canonical:
                    canonical = next((t for t in transcripts if t.get("translation")), None)
                if not canonical:
                    canonical = transcripts[0]

                translation = canonical.get("translation") or {}
                return {
                    "transcript_id": canonical["id"],
                    "chromosome": canonical.get("seq_region_name"),
                    "strand": canonical.get("strand"),
                    "exons": canonical.get("Exon", []),
                    "cds_start": translation.get("start"),
                    "cds_end": translation.get("end"),
                }
        except Exception as e:
            logger.error(f"❌ Transcript metadata failed: {e}")
            return None

    async def _build_cds_from_ucsc(self, metadata: Dict) -> Optional[str]:
        """Build the CDS by fetching each coding exon from UCSC."""
        chrom = f"chr{metadata['chromosome']}"
        strand = metadata["strand"]
        cds_start = metadata["cds_start"]
        cds_end = metadata["cds_end"]
        exons = metadata["exons"]

        # Sort coding exons in transcript order
        sorted_exons = sorted(exons, key=lambda e: e["start"], reverse=(strand == -1))

        fragments = []
        for exon in sorted_exons:
            e_start = max(exon["start"], cds_start)
            e_end = min(exon["end"], cds_end)
            if e_start >= e_end:
                continue
            # UCSC 0-based, half-open
            seq = await ucsc_service.fetch_region(chrom, e_start - 1, e_end)
            if not seq:
                return None
            fragments.append(seq)

        cds = "".join(fragments)
        if strand == -1:
            cds = ucsc_service.reverse_complement(cds)
        return cds

    async def _fetch_protein_from_ensembl(self, transcript_id: str) -> Optional[str]:
        """Fetch the protein sequence from Ensembl."""
        try:
            async with httpx.AsyncClient(timeout=30.0) as client:
                r = await client.get(
                    f"{self.base_url}/sequence/id/{transcript_id}?type=protein",
                    headers={"Content-Type": "application/json"},
                )
                if r.status_code == 200:
                    seq = r.json().get("seq")
                    if seq:
                        logger.info(f"✅ Protein from Ensembl: {len(seq)} aa")
                        return seq
        except Exception as e:
            logger.warning(f"⚠️ Ensembl protein failed: {e}")
        return None

    async def get_cds_sequence(self, gene: str) -> Optional[Dict]:
        if gene in self.cds_cache:
            return self.cds_cache[gene]

        # 1. Get metadata from Ensembl
        metadata = await self._fetch_transcript_metadata(gene)
        if not metadata:
            logger.warning(f"⚠️ Ensembl metadata failed for {gene}, using fallback")
            return self._get_fallback_cds(gene)

        # 2. Build CDS from UCSC
        cds = await self._build_cds_from_ucsc(metadata)
        if not cds:
            logger.warning(f"⚠️ UCSC CDS build failed for {gene}, using fallback")
            return self._get_fallback_cds(gene)

        # 3. Get protein from Ensembl
        protein = await self._fetch_protein_from_ensembl(metadata["transcript_id"])
        if not protein:
            logger.warning(f"⚠️ Ensembl protein failed, translating locally")
            protein = self.translate_cds(cds)

        result = {
            "cds": cds,
            "protein": protein,
            "transcript_id": metadata["transcript_id"],
            "gene": gene,
            "source": "ucsc_dna_ensembl_protein",
        }
        self.cds_cache[gene] = result
        logger.info(f"✅ CDS for {gene}: {len(cds)} bp, protein: {len(protein)} aa")
        return result

    def _get_fallback_cds(self, gene: str) -> Optional[Dict]:
        if gene in FALLBACK_CDS:
            logger.warning(f"⚠️ Using hardcoded fallback CDS for {gene}")
            result = FALLBACK_CDS[gene].copy()
            result["gene"] = gene
            result["source"] = "fallback"
            self.cds_cache[gene] = result
            return result
        return None

    def translate_cds(self, cds: str) -> str:
        # ... same as before ...
        ...


# Singleton
ensembl_service = EnsemblService()