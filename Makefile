lint:
	pre-commit run --all-files

publish:
	scripts/publish.sh

test:
	pytest -sv -n 4 tests/
	prove -v tests/*.t
