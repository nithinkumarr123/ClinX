"""ACMG/AMP 2015 classification."""

from enum import Enum
from typing import Dict, Tuple, Optional, List
from loguru import logger
import numpy as np


class ACMGClassification(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Variant of Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


class ACMGService:
    def __init__(self):
        self.criteria = {
            'PVS1': 1.0, 'PS1': 0.9, 'PS4': 0.9,
            'PM2': 0.6, 'PM3': 0.6, 'PP3': 0.3,
            'BS4': -0.9, 'BP5': -0.3, 'BP7': -0.3,
        }

    def classify(self, mutation_type, dna_delta: float, protein_delta: float,
                 clinvar_data: Optional[Dict] = None):
        evidence: Dict[str, List[str]] = {}
        total = 0.0
        mtype = mutation_type.value if hasattr(mutation_type, 'value') else str(mutation_type)

        if clinvar_data:
            sig = (clinvar_data.get('clinical_significance') or '').lower()
            if sig == 'benign':
                evidence['BS4'] = ["ClinVar: Benign"]
                return ACMGClassification.BENIGN, 0.05, evidence, "Benign (ClinVar)"
            if 'likely benign' in sig:
                evidence['BP5'] = ["ClinVar: Likely Benign"]
                return ACMGClassification.LIKELY_BENIGN, 0.15, evidence, "Likely Benign (ClinVar)"
            if sig == 'pathogenic':
                evidence['PS4'] = ["ClinVar: Pathogenic"]
                total += self.criteria['PS4']
            if 'likely pathogenic' in sig:
                evidence['PM3'] = ["ClinVar: Likely Pathogenic"]
                total += self.criteria['PM3']

        if mtype in ('frameshift', 'nonsense'):
            evidence['PVS1'] = ["Null variant - loss of function"]
            total += self.criteria['PVS1']

        if mtype == 'missense':
            if protein_delta < -5.0:
                evidence['PS1'] = [f"Severe protein impact (Δ={protein_delta:.3f})"]
                total += self.criteria['PS1']
            elif protein_delta < -2.0:
                evidence['PM2'] = [f"Moderate protein impact (Δ={protein_delta:.3f})"]
                total += self.criteria['PM2']
            elif protein_delta < -1.0:
                evidence['PP3'] = [f"Mild protein impact (Δ={protein_delta:.3f})"]
                total += self.criteria['PP3']

        if dna_delta < -0.5:
            evidence['PP3'] = [f"DNA disruptive (Δ={dna_delta:.3f})"]
            total += self.criteria['PP3']

        if mtype == 'synonymous':
            evidence['BP7'] = ["Synonymous variant"]
            total += self.criteria['BP7']

        posterior = 1.0 / (1.0 + np.exp(-total * 0.5))
        posterior = min(0.99, max(0.01, posterior))

        if total >= 1.0:
            cls, summ = ACMGClassification.PATHOGENIC, f"Pathogenic (posterior {posterior*100:.1f}%)"
        elif total >= 0.5:
            cls, summ = ACMGClassification.LIKELY_PATHOGENIC, f"Likely Pathogenic (posterior {posterior*100:.1f}%)"
        elif total <= -1.0:
            cls, summ = ACMGClassification.BENIGN, f"Benign (posterior {(1-posterior)*100:.1f}%)"
        elif total <= -0.5:
            cls, summ = ACMGClassification.LIKELY_BENIGN, f"Likely Benign (posterior {(1-posterior)*100:.1f}%)"
        else:
            cls, summ = ACMGClassification.VUS, f"VUS (posterior {posterior*100:.1f}%)"

        logger.info(f"📋 ACMG: {cls.value} (score={total:.2f})")
        return cls, posterior, evidence, summ


acmg_service = ACMGService()