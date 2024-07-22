import json
from pathlib import Path


# define function to read a json file

def read_json_file(file_path: Path):
    with open(file_path, 'r') as file:
        return json.load(file)

json_input_dir= Path('~/dev/madness_worktrees/test_madqc/json_input')


new_style_input= Path('~/dev/madness_worktrees/test_madqc/json_input_v2') 
# if new_style_input does not exist, create it

if not new_style_input.exists():
    new_style_input.mkdir()

#




