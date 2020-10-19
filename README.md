# scripts

Scripts organised by category.

Move this repo to `~/.local/scripts` and the following command to your shell rc file:
```
export PATH="$(find "$HOME/.local/scripts" -type d | grep -v "^.$\|.git\|pycache" | tr '\n' ':' | sed 's/:$//'):$PATH"
```
