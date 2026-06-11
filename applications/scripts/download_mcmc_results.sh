#!/bin/bash
# =============================================================================
# download_mcmc_results.sh
# Laedt Hu-Sawicki MCMC-Ergebnisse vom ellmos-services Server herunter.
#
# Server:  root@46.62.243.71 (ellmos-services)
# SSH-Key: ~/.ssh/id_ed25519_mcmc
# Quelle:  /opt/fst_calculations/husawicki_mcmc/
# Ziel:    .RESEARCH/Natur&Technik/3 Folgebeweise/Dark Energy/results/
#
# Erstellt: 2026-03-23
# =============================================================================

set -euo pipefail

# --- Konfiguration ---
SERVER="root@46.62.243.71"
SSH_KEY="$HOME/.ssh/id_ed25519_mcmc"
REMOTE_DIR="/opt/fst_calculations/husawicki_mcmc"
LOCAL_DIR="/c/Users/User/OneDrive/.RESEARCH/Natur&Technik/3 Folgebeweise/Dark Energy/results"

# Dateien die heruntergeladen werden sollen
FILES=(
    "husawicki_mcmc_chain.npz"
    "husawicki_mcmc.log"
    "mcmc_stdout.log"
    "compute_husawicki_mcmc.py"
)

# Beschreibungen fuer die Zusammenfassung
DESCRIPTIONS=(
    "MCMC-Chain (Samples, Log-Prob, Akzeptanzrate)"
    "MCMC-Logdatei (Parameter, Konvergenz)"
    "Stdout/Stderr des Rechenlaufs"
    "Python-Script (Backup der ausgefuehrten Version)"
)

# --- Farben ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# --- Funktionen ---
log_info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# --- Hauptprogramm ---
echo "============================================="
echo "  Hu-Sawicki MCMC -- Ergebnis-Download"
echo "============================================="
echo ""

# 1. Zielverzeichnis erstellen
if [ ! -d "$LOCAL_DIR" ]; then
    log_info "Erstelle Zielverzeichnis: $LOCAL_DIR"
    mkdir -p "$LOCAL_DIR"
fi

# 2. SSH-Verbindung testen
log_info "Teste SSH-Verbindung zu $SERVER ..."
if ! ssh -i "$SSH_KEY" -o ConnectTimeout=10 -o BatchMode=yes "$SERVER" "echo OK" > /dev/null 2>&1; then
    log_error "SSH-Verbindung fehlgeschlagen. Ist der Server erreichbar?"
    log_error "  ssh -i $SSH_KEY $SERVER"
    exit 1
fi
log_info "SSH-Verbindung OK."

# 3. Pruefen ob MCMC fertig ist (mcmc_stdout.log auf "completed" oder Chain-Datei vorhanden)
log_info "Pruefe ob MCMC abgeschlossen ist ..."
CHAIN_EXISTS=$(ssh -i "$SSH_KEY" "$SERVER" "test -f $REMOTE_DIR/husawicki_mcmc_chain.npz && echo yes || echo no")
if [ "$CHAIN_EXISTS" = "no" ]; then
    log_warn "Chain-Datei existiert noch nicht auf dem Server!"
    log_warn "MCMC laeuft moeglicherweise noch. Trotzdem verfuegbare Dateien herunterladen?"
    read -p "Fortfahren? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "Abgebrochen."
        exit 0
    fi
fi

# 4. Dateien herunterladen
echo ""
log_info "Starte Download ..."
echo ""

DOWNLOADED=0
FAILED=0
TOTAL_SIZE=0

for i in "${!FILES[@]}"; do
    FILE="${FILES[$i]}"
    DESC="${DESCRIPTIONS[$i]}"
    REMOTE_PATH="$REMOTE_DIR/$FILE"
    LOCAL_PATH="$LOCAL_DIR/$FILE"

    echo -n "  [$((i+1))/${#FILES[@]}] $FILE ... "

    # Pruefen ob Datei auf Server existiert
    FILE_EXISTS=$(ssh -i "$SSH_KEY" "$SERVER" "test -f $REMOTE_PATH && echo yes || echo no")
    if [ "$FILE_EXISTS" = "no" ]; then
        echo -e "${YELLOW}NICHT VORHANDEN (uebersprungen)${NC}"
        continue
    fi

    # Herunterladen
    if scp -i "$SSH_KEY" -q "$SERVER:$REMOTE_PATH" "$LOCAL_PATH" 2>&1; then
        # Dateigroesse pruefen
        if [ -f "$LOCAL_PATH" ]; then
            FILE_SIZE=$(stat --printf="%s" "$LOCAL_PATH" 2>/dev/null || stat -f%z "$LOCAL_PATH" 2>/dev/null || echo "0")
            if [ "$FILE_SIZE" -gt 0 ]; then
                # Menschenlesbare Groesse
                if [ "$FILE_SIZE" -gt 1048576 ]; then
                    SIZE_HR="$(( FILE_SIZE / 1048576 )) MB"
                elif [ "$FILE_SIZE" -gt 1024 ]; then
                    SIZE_HR="$(( FILE_SIZE / 1024 )) KB"
                else
                    SIZE_HR="$FILE_SIZE B"
                fi
                echo -e "${GREEN}OK${NC} ($SIZE_HR)"
                DOWNLOADED=$((DOWNLOADED + 1))
                TOTAL_SIZE=$((TOTAL_SIZE + FILE_SIZE))
            else
                echo -e "${RED}LEER (0 Bytes!)${NC}"
                FAILED=$((FAILED + 1))
            fi
        else
            echo -e "${RED}FEHLER (Datei nicht geschrieben)${NC}"
            FAILED=$((FAILED + 1))
        fi
    else
        echo -e "${RED}FEHLER (scp fehlgeschlagen)${NC}"
        FAILED=$((FAILED + 1))
    fi
done

# 5. Zusammenfassung
echo ""
echo "============================================="
echo "  Zusammenfassung"
echo "============================================="
echo ""

if [ "$TOTAL_SIZE" -gt 1048576 ]; then
    TOTAL_HR="$(( TOTAL_SIZE / 1048576 )) MB"
elif [ "$TOTAL_SIZE" -gt 1024 ]; then
    TOTAL_HR="$(( TOTAL_SIZE / 1024 )) KB"
else
    TOTAL_HR="$TOTAL_SIZE B"
fi

log_info "Heruntergeladen:  $DOWNLOADED / ${#FILES[@]} Dateien"
log_info "Gesamtgroesse:    $TOTAL_HR"
log_info "Zielverzeichnis:  $LOCAL_DIR"

if [ "$FAILED" -gt 0 ]; then
    log_warn "Fehlgeschlagen:   $FAILED Dateien"
fi

echo ""

# 6. Inhalt des Zielverzeichnisses anzeigen
log_info "Inhalt des Zielverzeichnisses:"
echo ""
ls -lh "$LOCAL_DIR/" 2>/dev/null || log_warn "Verzeichnis leer oder nicht lesbar."
echo ""

# 7. Kurzanalyse der Chain-Datei (falls vorhanden und Python verfuegbar)
CHAIN_FILE="$LOCAL_DIR/husawicki_mcmc_chain.npz"
if [ -f "$CHAIN_FILE" ]; then
    log_info "Kurzanalyse der MCMC-Chain:"
    echo ""
    PYTHONIOENCODING=utf-8 python -c "
import numpy as np
import sys

try:
    data = np.load('${CHAIN_FILE//\'/\\\'}')
    print(f'  Enthaltene Arrays: {list(data.keys())}')
    for key in data.keys():
        arr = data[key]
        print(f'  {key}: shape={arr.shape}, dtype={arr.dtype}')
        if arr.ndim >= 1 and arr.shape[0] > 0:
            if np.issubdtype(arr.dtype, np.floating):
                print(f'    min={arr.min():.6f}, max={arr.max():.6f}, mean={arr.mean():.6f}')
    print()
    print('  Chain-Datei erfolgreich gelesen.')
except Exception as e:
    print(f'  Fehler beim Lesen: {e}', file=sys.stderr)
" 2>&1 || log_warn "Python-Analyse fehlgeschlagen (numpy installiert?)."
    echo ""
fi

# 8. Naechste Schritte
echo "============================================="
echo "  Naechste Schritte"
echo "============================================="
echo ""
echo "  1. Chain analysieren:  python -c \"import numpy as np; d=np.load('$CHAIN_FILE'); ...\""
echo "  2. Konvergenz pruefen: Trace-Plots, Gelman-Rubin R-hat, effektive Samplegroesse"
echo "  3. Ergebnisse in Paper einfuegen: FST-DE_DarkEnergy_Skeleton_v1_{en,de}.tex"
echo "  4. Bei Bedarf Server aufraeumen: ssh -i $SSH_KEY $SERVER 'rm -rf $REMOTE_DIR'"
echo ""

if [ "$FAILED" -eq 0 ] && [ "$DOWNLOADED" -gt 0 ]; then
    log_info "Download abgeschlossen. Alle Dateien OK."
    exit 0
elif [ "$DOWNLOADED" -eq 0 ]; then
    log_error "Keine Dateien heruntergeladen!"
    exit 1
else
    log_warn "Download mit Warnungen abgeschlossen ($FAILED Fehler)."
    exit 1
fi
