import glob
import re

def replace_tone(text):
    # DEPRECATED: This function was the source of manuscript inflation by
    # automatically replacing honest scientific language with bullish prose.
    # It is now a no-op to preserve accurate terminology.
    #
    # Previously it replaced:
    #   "noise-limited" → "ambient screening verified"
    #   "diagnostic" → "validation"
    #   etc.
    #
    # These replacements overclaimed non-detections as confirmations and
    # misrepresented the epistemic status of supporting evidence.
    return text

def main():
    print("refactor_html.py is DEPRECATED and no longer modifies HTML files.")
    print("Honest scientific terminology is preserved in site/components/.")

if __name__ == "__main__":
    main()
