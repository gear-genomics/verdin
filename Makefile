TARGETS = .silica

all:   	$(TARGETS)
	cd src/silica && make all && cd ../../ && touch .silica

clean:
	cd src/silica && make clean
	rm -f $(TARGETS) $(TARGETS:=.o)
