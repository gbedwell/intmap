import json
from pathlib import Path

def load_test_data(module, test_case):
    data_dir = Path(__file__).parent / 'data' / module / test_case
    
    with open(data_dir / 'input.json') as f:
        input_data = json.load(f)
    with open(data_dir / 'expected.json') as f:
        expected = json.load(f)
        
    return input_data, expected