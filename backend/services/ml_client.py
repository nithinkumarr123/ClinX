"""ML Client - Communicates with Colab ML server."""

import httpx
from typing import Dict
from loguru import logger
from config import config


class MLClient:
    def __init__(self):
        self.base_url = (config.COLAB_API_URL or "").rstrip("/")
        self.timeout = config.API_TIMEOUT
        self.headers = {
            "Content-Type": "application/json",
            "ngrok-skip-browser-warning": "true",
            "User-Agent": "ClinX-Backend/1.0",
        }
        logger.info(f"ML Client initialized with base URL: {self.base_url}")

    async def score_dna(self, ref_seq: str, alt_seq: str, position: int) -> Dict:
        if not self.base_url:
            return {"delta_score": 0.0, "interpretation": "not_configured"}
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                r = await client.post(
                    f"{self.base_url}/score_dna",
                    json={"ref": ref_seq, "alt": alt_seq, "pos": position},
                    headers=self.headers,
                )
                r.raise_for_status()
                return r.json()
        except Exception as e:
            logger.error(f"DNA scoring failed: {e}")
            return {"delta_score": 0.0, "interpretation": "error"}

    async def score_protein(self, sequence: str, position: int, ref_aa: str, alt_aa: str) -> Dict:
        if not self.base_url:
            return {"delta_score": 0.0, "impact": "not_configured"}
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                r = await client.post(
                    f"{self.base_url}/score_protein",
                    json={"sequence": sequence, "position": position, "ref_aa": ref_aa, "alt_aa": alt_aa},
                    headers=self.headers,
                )
                r.raise_for_status()
                return r.json()
        except Exception as e:
            logger.error(f"Protein scoring failed: {e}")
            return {"delta_score": 0.0, "impact": "error"}

    async def health_check(self) -> Dict:
        if not self.base_url:
            return {"status": "not_configured"}
        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                r = await client.get(f"{self.base_url}/health", headers=self.headers)
                return r.json() if r.status_code == 200 else {"status": "unhealthy"}
        except Exception as e:
            return {"status": "unreachable", "error": str(e)}


ml_client = MLClient()