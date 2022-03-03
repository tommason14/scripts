#!/usr/bin/env bash

# read command is different in bash and zsh - using bash form here
read -p "Enter language or core util: " first
read -p "Enter topic (press enter if core util): " second

if [[ $second ]] ; then
  curl -s cht.sh/$first/${second// /+} | less -R
else
  curl -s cht.sh/$first | less -R
fi
