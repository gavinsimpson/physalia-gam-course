.PHONY : slides
slides:
	cd ./day-1 && make $@
	cd ./day-2 && make $@
	cd ./day-3 && make $@
	cd ./day-4 && make $@
	cd ./day-5 && make $@

.PHONY : purl
purl:
	cd ./dat-1 && make $@
	cd ./dat-2 && make $@
	cd ./dat-3 && make $@
	cd ./day-4 && make $@
	cd ./day-5 && make $@

