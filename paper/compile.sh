woodchuck u_vel.csv u_vel.tex \
	--caption "Velocity in the x direction of the solver for two different reynolds numbers compared to truth values from literature" \
	--label "table_u_vel" \
	--table-args "" \
	--table-format "*"
woodchuck v_vel.csv v_vel.tex \
	--caption "Velocity in the y direction of the solver for two different reynolds numbers compared to truth values from literature" \
	--label "table_v_vel" \
	--table-args "" \
	--table-format "*"
pdflatex report && bibtex report && pdflatex report
