.PHONY: build run up down clean

build:
	docker compose build

run:
	docker compose up
	Rscript src/main.R
	Rscript src/update_readme.R

up:
	docker compose up --build

down:
	docker compose down

clean:
	docker compose down --rmi local
	find data/cln -type f -delete
