#!/bin/bash
# GitHub Setup Script for DGAT_Immunology
# 
# BEFORE RUNNING THIS:
# 1. Create repository on GitHub: https://github.com/new
#    - Name: DGAT_Immunology
#    - Description: Multi-database analysis of DGAT1 expression and immune microenvironment in glioblastoma
#    - Public (recommended) or Private
#    - DO NOT initialize with README, .gitignore, or license
# 
# 2. Copy your new repository URL (it will look like):
#    git@github.com:ZijieFeng-98/DGAT_Immunology.git
# 
# 3. Replace YOUR_REPO_URL below with your actual URL

# =====================================================
# CONFIGURATION - EDIT THIS!
# =====================================================
YOUR_REPO_URL="git@github.com:ZijieFeng-98/DGAT_Immunology.git"
# =====================================================

cd "/Users/zijiefeng/Desktop/Guo's lab/My_Research/DGAT_Immunology"

echo "🚀 Setting up GitHub repository for DGAT_Immunology..."
echo ""

# Add remote
echo "📡 Adding remote repository..."
git remote add origin "$YOUR_REPO_URL"
echo "✅ Remote added"
echo ""

# Verify remote
echo "🔍 Verifying remote..."
git remote -v
echo ""

# Stage new files
echo "📦 Staging files..."
git add .gitignore
git add Notes/Daily_Logs/
git add DATA_GAP_ANALYSIS.md
git add DATA_STATUS_SUMMARY.md
git add DATA_CHECK_RESULTS.txt
git add PROJECT_STANDARDS_CHECKLIST.md
git add STANDARDS_BASELINE_SUMMARY.md
git add GITHUB_SETUP_NEW_REPO.md
git add Processed_Data/*/README.md
git add Processed_Data/*/*.txt
git add -u  # Stage deletions
echo "✅ Files staged"
echo ""

# Show what will be committed
echo "📋 Files to be committed:"
git status --short
echo ""

# Commit
echo "💾 Creating commit..."
git commit -m "Project setup: standards baseline, data analysis pipeline, documentation

- Add comprehensive bioinformatics standards checklist (Nature/Cell requirements)
- Add data gap analysis and integration plan
- Add daily logging system
- Clean up unnecessary README files
- Update .gitignore to exclude large data files
- Document Week 1 priorities and current status

Status: 70% publication-ready, metadata integration needed"
echo "✅ Commit created"
echo ""

# Ensure main branch
echo "🌿 Setting main branch..."
git branch -M main
echo "✅ Main branch set"
echo ""

# Push
echo "🚀 Pushing to GitHub..."
echo "   This may take a few moments..."
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ ✅ ✅ SUCCESS! ✅ ✅ ✅"
    echo ""
    echo "Your repository is now live at:"
    echo "   https://github.com/$(echo $YOUR_REPO_URL | sed 's/.*github.com[:/]//' | sed 's/.git$//')"
    echo ""
    echo "Next steps:"
    echo "  1. Visit your repository on GitHub"
    echo "  2. Add topics: glioblastoma, cancer-immunology, tcga, bioinformatics"
    echo "  3. Consider adding a LICENSE (MIT recommended)"
    echo ""
else
    echo ""
    echo "❌ Push failed. Common issues:"
    echo "   1. Make sure you created the repository on GitHub first"
    echo "   2. Check that YOUR_REPO_URL is correct in this script"
    echo "   3. Verify SSH key is set up: ssh -T git@github.com"
    echo ""
    echo "See GITHUB_SETUP_NEW_REPO.md for troubleshooting"
fi

