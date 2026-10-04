# Überholte Ergebnisse

`SHELL_DFC_WATERLINE_LEDGER_2026-06-05.{csv,json,md}` wurden mit der fehlerhaften
Sabra-Dynamik (vor Issue #1) und der dazu gehörenden Flussformel berechnet. Die Dateien
bleiben zur Nachvollziehbarkeit erhalten, sind aber **überholt**; ihre Kennzahlen
(z. B. `weighted/target=0.217790`, `residual/allowance=5.21473`) dürfen nicht
mehr zitiert werden.

Nachfolger (korrigierte Dynamik, energiekonsistenter Fluss `Pi_n`):
`SHELL_DFC_WATERLINE_LEDGER_2026-10-04_sabra-v2.{csv,json,md}`, erzeugt von
`../compute_shell_dfc_waterline_ledger.py` (`MODEL_TAG = "sabra-v2"`).

`DUAL_DFC_LEDGER_2026-05-27.*` ist ein Toy-Ledger ohne Shell-Dynamik und davon nicht betroffen.

Zahlen im Paper (`fst-physics/turbulence/FST-TU_Turbulence_Skeleton_v1_en.tex`) müssen
mit der korrigierten Dynamik neu validiert werden.
