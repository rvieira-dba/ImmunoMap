"""Allow ``python -m immunomap``."""
import sys

from immunomap.cli.main import main

if __name__ == "__main__":
    sys.exit(main())
