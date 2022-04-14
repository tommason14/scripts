#!/usr/bin/env bash

read -p "Search term: " search_term
bibsearch search "$search_term" | less 
