.PHONY: install test run lint
install:
	python -m pip install -e . -r requirements.txt
test:
	python -m pytest
run:
	python examples/argon_trimer.py
lint:
	python -m compileall -q src examples
