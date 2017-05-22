" Vim syntax file 
" Language : Bathe Input File
" Maintainer : Guido Giuntoli

au BufRead,BufNewFile *.bathe set filetype=bathe

if exists("b:current_syntax")
      finish
endif

" Keywords
syn keyword syntaxElementKeyword  $Mesh $EndMesh $Materials $EndMaterials  nextgroup=syntaxElement2

" Matches
syn match syntaxElementMatch 'regexp' contains=syntaxElement1 nextgroup=syntaxElement2 skipwhite

" Regions
syn region syntaxElementRegion start='x' end='y'
