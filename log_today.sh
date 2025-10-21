#!/bin/bash
# =============================================================================
# log_today.sh - Quick logging wrapper
# =============================================================================
#
# Usage:
#   ./log_today.sh "Started analyzing TCGA data"
#   ./log_today.sh --view
#   ./log_today.sh --weekly
# =============================================================================

# Navigate to project directory
cd "$(dirname "$0")"

if [ "$1" = "--view" ]; then
    Rscript Scripts/log_activity.R --view
elif [ "$1" = "--weekly" ]; then
    Rscript Scripts/log_activity.R --weekly
elif [ "$1" = "--summary" ]; then
    Rscript Scripts/log_activity.R --summary
elif [ -n "$1" ]; then
    # Log a message
    Rscript -e "source('Scripts/log_activity.R'); log_activity('$*', category='general')"
else
    # Open today's log in default editor
    LOG_FILE="Notes/Daily_Logs/WorkLog_$(date +%Y-%m-%d).md"
    if [ -f "$LOG_FILE" ]; then
        open "$LOG_FILE"
    else
        echo "No log file for today. Creating..."
        Rscript -e "source('Scripts/log_activity.R'); log_activity('Started work session', category='general')"
        open "$LOG_FILE"
    fi
fi

