.PHONY: setup pipeline dashboard clean

PYTHON ?= python3

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) dashboard.py

clean:
	rm -f loblaw.db
	rm -rf outputs
	mkdir -p outputs