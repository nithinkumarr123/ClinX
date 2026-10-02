"""Calibration Service."""

from typing import Dict


class CalibrationService:
    def calibrate_dna_score(self, delta_score: float, gene: str) -> Dict:
        if delta_score < -0.5:
            interp = "disruptive"
        elif delta_score < -0.2:
            interp = "moderate"
        elif delta_score < 0.2:
            interp = "neutral"
        elif delta_score < 0.5:
            interp = "mildly tolerated"
        else:
            interp = "tolerated"
        return {"calibrated_delta": round(delta_score, 3), "interpretation": interp}


calibration_service = CalibrationService()