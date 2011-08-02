	: LOOP
	s/\\\$//g
/^[^$]*\$[^$]*$/ N
/\\cite{[^}]*,$/ N
/\\cite{[^}]*,$/ N
	s/\n/ /
	s/\\\[/\$/g
	s/\\\]/\$/g
	s/\\begin{eqnarray\*}/\$/g
	s/\\end{eqnarray\*}/\$/g

	s/\\begin{[^}]*}//g
	s/\\end{[^}]*}//g
	s/\\ref{[^}]*}//g
	s/\\cite{[^}]*}//g
	s/\\label{[^}]*}//g
	s/\\usepackage{[^}]*}//g
	s/\\bibliography{[^}]*}//g
	s/\\def\\[a-zA-Z]*{[^}]*}//g
	s/[a-zA-Z.]*@[a-zA-Z.]*//g
	s/Erd\\H os//g
	s/Erd\\H{o}s//g
	s/\\'//g

	s/\\[a-zA-Z]*/ /g

        s/\%.*$//
	s/\$[^$]\+\$//g
/^[^$]*\$[^$]*$/ t LOOP
