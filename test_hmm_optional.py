#!/usr/bin/env python3
"""Test script to verify optional HMM functionality works correctly."""

import sys
import tempfile
import os

def test_optional_hmm_import():
    """Test that HMM import is optional and provides helpful error messages."""
    # Temporarily hide pomegranate from imports
    original_modules = sys.modules.copy()
    
    # Remove pomegranate modules if they exist
    to_remove = [mod for mod in sys.modules if 'pomegranate' in mod or 'torch' in mod]
    for mod in to_remove:
        del sys.modules[mod]
    
    # Add a mock module that raises ImportError
    class MockPomegranate:
        def __getattr__(self, name):
            raise ImportError("No module named 'pomegranate'")
    
    sys.modules['pomegranate'] = MockPomegranate()
    
    try:
        # This should work without error, but HMM methods should be unavailable
        from cnvlib.segmentation import SEGMENT_METHODS, HMM_METHODS
        
        print(f"Available methods: {SEGMENT_METHODS}")
        print(f"HMM methods: {HMM_METHODS}")
        
        # HMM methods should be empty
        assert len(HMM_METHODS) == 0, f"Expected no HMM methods, got {HMM_METHODS}"
        
        # Basic methods should still be available
        expected_basic = ("cbs", "flasso", "haar", "none")
        for method in expected_basic:
            assert method in SEGMENT_METHODS, f"Expected {method} in {SEGMENT_METHODS}"
        
        print("✓ Optional HMM import works correctly")
        
        # Test error message when trying to use HMM
        try:
            from cnvlib.segmentation import do_segmentation
            from cnvlib.cnary import CopyNumArray as CNA
            import pandas as pd
            
            # Create a minimal test array
            test_data = pd.DataFrame({
                'chromosome': ['1'] * 3,
                'start': [1000, 2000, 3000],
                'end': [2000, 3000, 4000],
                'log2': [0.0, 0.5, -0.5],
                'weight': [1.0, 1.0, 1.0],
                'gene': ['A', 'B', 'C']
            })
            cnarr = CNA(test_data)
            
            # This should raise ImportError with helpful message
            do_segmentation(cnarr, "hmm")
            assert False, "Expected ImportError for HMM method"
            
        except ImportError as e:
            error_msg = str(e)
            assert "pomegranate >= 1.0.0" in error_msg, f"Error message should mention pomegranate: {error_msg}"
            assert "pip install cnvkit[hmm]" in error_msg, f"Error message should mention installation: {error_msg}"
            print("✓ Helpful error message provided for missing HMM dependencies")
            
    finally:
        # Restore original modules
        sys.modules.clear()
        sys.modules.update(original_modules)

if __name__ == "__main__":
    test_optional_hmm_import()
    print("All tests passed!")