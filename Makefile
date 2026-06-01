.PHONY: build run up down clean

build:
	docker compose build

run:
	docker compose up

up:
	docker compose up --build

down:
	docker compose down

clean:
	docker compose down --rmi local
	find data/cln -type f -delete
