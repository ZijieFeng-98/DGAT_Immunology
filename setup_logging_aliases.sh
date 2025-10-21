#!/bin/bash
# =============================================================================
# setup_logging_aliases.sh - Add convenient logging aliases to your shell
# =============================================================================
#
# This adds shortcuts to log your work from anywhere in the terminal
#
# Usage: bash setup_logging_aliases.sh
# =============================================================================

PROJECT_DIR="/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"
SHELL_RC="$HOME/.zshrc"

echo "Setting up daily logging aliases..."
echo ""

# Create alias definitions
ALIASES="
# ============================================
# Daily Work Log Aliases (Added $(date))
# ============================================
alias logwork='cd \"$PROJECT_DIR\" && ./log_today.sh'
alias viewlog='cd \"$PROJECT_DIR\" && ./log_today.sh --view'
alias weeklog='cd \"$PROJECT_DIR\" && ./log_today.sh --weekly'
alias dgat='cd \"$PROJECT_DIR\"'
"

# Check if aliases already exist
if grep -q "Daily Work Log Aliases" "$SHELL_RC" 2>/dev/null; then
    echo "⚠️  Aliases already exist in $SHELL_RC"
    echo "   Skipping to avoid duplicates."
else
    # Add aliases to shell config
    echo "$ALIASES" >> "$SHELL_RC"
    echo "✅ Added aliases to $SHELL_RC"
fi

echo ""
echo "📝 Aliases added:"
echo "   logwork 'your message'  - Log from anywhere"
echo "   viewlog                 - View today's log"
echo "   weeklog                 - View weekly summary"
echo "   dgat                    - Jump to project directory"
echo ""
echo "🔄 To activate now, run:"
echo "   source ~/.zshrc"
echo ""
echo "Then try:"
echo "   logwork 'Set up automated logging!'"
echo ""

