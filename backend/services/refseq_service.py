"""
RefSeq Service - Fetch CDS and protein directly from NCBI.
This avoids all exon-assembly issues.
"""

import httpx
from typing import Optional, Dict
from loguru import logger
from tenacity import retry, stop_after_attempt, wait_exponential
from config import config


class RefSeqService:
    def __init__(self):
        self.base_url = config.NCBI_API
        self.cache: Dict[str, Dict] = {}

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=5))
    async def _fetch_fasta(self, db: str, accession: str, rettype: str) -> Optional[str]:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(
                f"{self.base_url}/efetch.fcgi",
                params={
                    "db": db,
                    "id": accession,
                    "rettype": rettype,
                    "retmode": "text",
                },
            )
            r.raise_for_status()
            lines = r.text.strip().split("\n")
            seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            return seq.upper() if seq else None

    async def get_cds_and_protein(self, gene: str, refseq_mrna: str, refseq_protein: str) -> Optional[Dict]:
        if gene in self.cache:
            return self.cache[gene]

        logger.info(f"🌐 Fetching {refseq_mrna} CDS from NCBI...")
        # fasta_cds_na returns the CDS directly (not the mRNA)
        cds = await self._fetch_fasta("nuccore", refseq_mrna, "fasta_cds_na")
        if not cds:
            logger.error(f"❌ Could not fetch CDS for {refseq_mrna}")
            return None

        logger.info(f"🌐 Fetching {refseq_protein} from NCBI...")
        protein = await self._fetch_fasta("protein", refseq_protein, "fasta")

        result = {
            "cds": cds,
            "protein": protein,
            "transcript_id": refseq_mrna,
            "protein_id": refseq_protein,
            "gene": gene,
            "source": "ncbi_refseq",
        }
        self.cache[gene] = result
        logger.info(f"✅ {gene}: CDS={len(cds)} bp, protein={len(protein) if protein else 0} aa")
        return result

    async def health(self) -> bool:
        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                r = await client.get(f"{self.base_url}/einfo.fcgi", params={"retmode": "json"})
                return r.status_code == 200
        except Exception:
            return False


refseq_service = RefSeqService()