from pathlib import Path

from pypdf import PdfWriter


ROOT = Path(__file__).resolve().parent
EN_PDF = ROOT / "PvsNP_Entropy_EN.pdf"
DE_PDF = ROOT / "PvsNP_Entropy_DE.pdf"
KOMBI_PDF = ROOT / "PvsNP_Entropy_kombi.pdf"


def main() -> None:
    writer = PdfWriter()
    writer.append(str(EN_PDF), import_outline=False)
    writer.append(str(DE_PDF), import_outline=False)
    writer.add_metadata(
        {
            "/Title": "An Entropic Perspective on P vs NP / Eine entropische Perspektive auf P vs NP",
            "/Author": "Lukas Geiger",
            "/Subject": "Combined bilingual edition of the P vs NP entropy paper",
            "/Keywords": (
                "P vs NP, Kolmogorov complexity, algorithmic entropy, witness compression, "
                "bilingual combined edition"
            ),
        }
    )
    with KOMBI_PDF.open("wb") as fh:
        writer.write(fh)


if __name__ == "__main__":
    main()
