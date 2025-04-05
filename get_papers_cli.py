
import click
from get_papers.pubmed_client import fetch_papers
import csv

@click.command()
@click.argument("query")
@click.option("-d", "--debug", is_flag=True, help="Enable debug mode")
@click.option("-f", "--file", type=click.Path(), help="Output CSV file")
def main(query: str, debug: bool, file: str | None):
    papers = fetch_papers(query, debug=debug)

    if file:
        with open(file, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=papers[0].keys())
            writer.writeheader()
            writer.writerows(papers)
        click.echo(f"Saved {len(papers)} results to {file}")
    else:
        for paper in papers:
            click.echo(paper)

if __name__ == "__main__":
    main()
