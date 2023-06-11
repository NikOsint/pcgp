.PHONY: run
run: pcgp.exe
	@if [ -z "$(N)" ]; then echo "Parameter N not difined"; \
	else \
		if [ -z "$(K)" ]; then echo "Parameter K not difined"; \
		else \
			if [ -z "$(C)" ]; then echo "Parameter C not difined"; \
			else \
				rm -rf pcgp-$(N)-$(K)-$(C) && mkdir pcgp-$(N)-$(K)-$(C) && cd pcgp-$(N)-$(K)-$(C) && ../pcgp.exe s $(N) $(K) $(C); \
			fi \
		fi \
	fi

pcgp.exe: pcgpmem.h pcgpmem.cpp pcgp.cpp
	g++ -Wall -std=c++11 pcgpmem.cpp pcgp.cpp -o pcgp.exe

.PHONY: clean
clean:
	rm -rf pcgp.exe pcgp-*/

