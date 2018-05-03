STRUCT_SRC=$(wildcard structures/*.sdf)
LIGPREP_FILES=$(patsubst structures/%.sdf, ligprep/%.sdf, $(STRUCT_SRC))
PROP_FILES=$(patsubst structures/%.sdf, properties/%.csv, $(STRUCT_SRC))
SCHRODINGER=/opt/schrodinger/suites2018-1
LIGPREP_BIN=$(SCHRODINGER)/ligprep
LIGPREP_OPTS=-inp ../../ligprep.inp -WAIT -NJOBS 4 -HOST localhost:2 -LOCAL

.PHONY : all wash_molecules calc_props
all : calc_props wash_molecules

wash_molecules : $(LIGPREP_FILES)
ligprep/%.sdf : structures/%.sdf
	mkdir -p ligprep/$*-job
	cd ligprep/$*-job && $(LIGPREP_BIN) $(LIGPREP_OPTS) -isd ../../$< -osd ../../$@ -TMPDIR .

calc_props : $(PROP_FILES)
properties/%.csv : ligprep/%.sdf calc_props.py
	python calc_props.py -i $< -o properties -c

test:
	nosetests tests

clean :
	rm properties/*.csv ligprep/*.sdf *.log
