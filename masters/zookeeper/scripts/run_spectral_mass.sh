#!/bin/bash
# Wartet bis der aktuelle Lauf (PID 1399673) fertig ist, dann startet spectral_mass
echo "Warte auf PID 1399673 (c2bt_server_quasimode_scaling.py)..."
while kill -0 1399673 2>/dev/null; do
    sleep 30
done
echo "Erster Lauf fertig. Starte Spektralmassen-Analyse..."
cd /root/c2bt_run
PYTHONIOENCODING=utf-8 DPS=35 N_FACTOR=12 python3 -u c2bt_spectral_mass.py > c2bt_spectral_mass_output.log 2>&1
echo "Spektralmassen-Analyse fertig."
# Ollama wieder starten
systemctl start ollama 2>/dev/null
echo "Ollama gestartet."
