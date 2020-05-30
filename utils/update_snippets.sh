#!/usr/bin/env bash

if [[ $USER =~ (tommason|tmas0023) ]]
then
  snipdir="$HOME/Documents/repos/vim-snippets"
else
  snipdir="$HOME/vim-snippets"
fi

# Push to github
cd "$snipdir"
git add .
git commit -m "Snippets updated"
git push 

# Pull changes into $HOME/.vim/bundle
vim -c ":PluginUpdate vim-snippets" -c ":q" -c ":q"
