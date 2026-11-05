"""Rebuild the integrated data release offline from immutable downloaded sources."""
from src.unify_datasets import run
from src.package_unified_dataset import write_report
if __name__=='__main__':
    run()
    write_report()
