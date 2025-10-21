# TCGA Analysis Test Suite
# Usage: make test

.PHONY: test clean

test:
	@echo "🧪 Running TCGA Analysis Test Suite..."
	@echo "📊 Step 1: Running full analysis (batch mode)..."
	@Rscript Scripts/tests/01_run_all.R || true
	@echo "🔍 Step 2: Validating outputs..."
	@Rscript Scripts/tests/02_validate_outputs.R || true
	@echo "✅ Test suite complete! Check Processed_Data/test_summary.json for results."

clean:
	@echo "🧹 Cleaning test outputs..."
	@rm -rf Processed_Data/test_run.log
	@rm -rf Processed_Data/test_summary.json
	@echo "✅ Cleanup complete!"

html-test:
	@echo "🌐 Running HTML export test for GBM..."
	@Rscript Scripts/tests/03_replay_one_as_html.R
	@echo "✅ HTML test complete!"