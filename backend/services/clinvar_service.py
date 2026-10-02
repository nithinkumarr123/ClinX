"""ClinVar Service."""

import httpx
from typing import Dict
from loguru import logger
from config import config


class ClinVarService:
    def __init__(self):
        self.base_url = config.NCBI_API
        self.known_variants = {
            ("TP53", 7675088, "G", "A"): {"sig": "Pathogenic", "cond": ["Li-Fraumeni syndrome"], "rev": "Reviewed by expert panel"},
            ("TP53", 7676152, "C", "G"): {"sig": "Benign", "cond": ["Not specified"], "rev": "Reviewed by expert panel"},
            ("TP53", 7674924, "G", "A"): {"sig": "Pathogenic", "cond": ["Li-Fraumeni syndrome"], "rev": "Reviewed by expert panel"},
            ("TP53", 7674240, "G", "A"): {"sig": "Pathogenic", "cond": ["Li-Fraumeni syndrome"], "rev": "Reviewed by expert panel"},
            ("TP53", 7674182, "C", "T"): {"sig": "Pathogenic", "cond": ["Li-Fraumeni syndrome"], "rev": "Reviewed by expert panel"},
            ("TP53", 7673632, "G", "A"): {"sig": "Pathogenic", "cond": ["Li-Fraumeni syndrome"], "rev": "Reviewed by expert panel"},
            ("BRCA1", 43091971, "AG", ""): {"sig": "Pathogenic", "cond": ["Hereditary breast ovarian cancer"], "rev": "Reviewed by expert panel"},
            ("BRCA1", 43091994, "C", "A"): {"sig": "Benign", "cond": ["Not specified"], "rev": "Reviewed by expert panel"},
            ("BRCA1", 43091994, "C", "G"): {"sig": "Benign", "cond": ["Not specified"], "rev": "Reviewed by expert panel"},
            ("BRCA1", 43098602, "C", "T"): {"sig": "Pathogenic", "cond": ["Hereditary breast ovarian cancer"], "rev": "Reviewed by expert panel"},
            ("BRCA1", 43090451, "T", "G"): {"sig": "Likely Pathogenic", "cond": ["Hereditary breast ovarian cancer"], "rev": "Reviewed by expert panel"},
        }

    async def get_clinvar_for_variant(self, gene: str, ref: str, alt: str, position: int = None) -> Dict:
        if position is not None:
            key = (gene, position, ref, alt)
            if key in self.known_variants:
                kv = self.known_variants[key]
                logger.info(f"✅ Known ClinVar: {gene}:{position} → {kv['sig']}")
                return {
                    "clinical_significance": kv["sig"],
                    "conditions": kv["cond"],
                    "review_status": kv["rev"],
                    "found": True,
                }

        try:
            async with httpx.AsyncClient(timeout=20.0) as client:
                term = f"{gene}[gene]"
                if position:
                    term += f" AND {position}[chrpos]"
                r = await client.get(
                    f"{self.base_url}/esearch.fcgi",
                    params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 5},
                )
                ids = r.json().get("esearchresult", {}).get("idlist", [])
                if ids:
                    r2 = await client.get(
                        f"{self.base_url}/esummary.fcgi",
                        params={"db": "clinvar", "id": ",".join(ids[:3]), "retmode": "json"},
                    )
                    result = r2.json().get("result", {})
                    sigs = []
                    for uid in ids[:3]:
                        doc = result.get(uid, {})
                        cs = doc.get("clinical_significance", "")
                        if isinstance(cs, dict):
                            cs = cs.get("description", "")
                        if cs:
                            sigs.append(cs)
                    sig = "Uncertain significance"
                    for cand in ["Pathogenic", "Likely pathogenic", "Benign", "Likely benign"]:
                        if any(cand.lower() in s.lower() for s in sigs):
                            sig = cand
                            break
                    return {
                        "clinical_significance": sig,
                        "conditions": ["Not specified"],
                        "review_status": "Available",
                        "found": True,
                    }
        except Exception as e:
            logger.warning(f"⚠️ ClinVar API failed: {e}")

        return {
            "clinical_significance": "Uncertain significance",
            "conditions": ["Not specified"],
            "review_status": "No review",
            "found": False,
        }


clinvar_service = ClinVarService()