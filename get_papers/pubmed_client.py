import ssl
import urllib.request
ssl._create_default_https_context = ssl._create_unverified_context

from typing import List, Dict
from Bio import Entrez
from xml.etree import ElementTree as ET

Entrez.email = "pratiksharkutal@gmail.com"  

def fetch_papers(query: str, debug: bool = False, max_results: int = 10) -> List[Dict]:
    if debug:
        print(f"[DEBUG] Querying PubMed for: {query}")

    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_results.get("IdList", [])
    if debug:
        print(f"[DEBUG] Found {len(id_list)} articles")

    if not id_list:
        return []

    fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
    data = fetch_handle.read()
    fetch_handle.close()

    root = ET.fromstring(data)
    papers = []

    for article in root.findall(".//PubmedArticle"):
        medline = article.find("MedlineCitation")
        article_data = medline.find("Article")

        title = article_data.findtext("ArticleTitle", default="No Title")
        pub_date = article_data.find(".//PubDate")

        pub_date_text = "Unknown"
        if pub_date is not None:
            year = pub_date.findtext("Year")
            month = pub_date.findtext("Month")
            day = pub_date.findtext("Day")
            pub_date_text = f"{year or ''}-{month or ''}-{day or ''}".strip("-")

        authors = article_data.findall("AuthorList/Author")
        non_academic_authors = []
        company_affiliations = []
        corresponding_email = None

        for author in authors:
            affiliation = author.findtext("AffiliationInfo/Affiliation")
            email = None
            if affiliation:
                
                affil_lower = affiliation.lower()
                if not any(word in affil_lower for word in ["university", "college", "institute", "school", "hospital", "clinic"]):
                    non_academic_authors.append(author.findtext("LastName", "Unknown"))
                    company_affiliations.append(affiliation)
                
                if "@" in affiliation and not corresponding_email:
                    email = affiliation.split()[-1]
                    if "@" in email:
                        corresponding_email = email.strip("().;,:")

        papers.append({
            "PubmedID": medline.attrib.get("PMID", "Unknown"),
            "Title": title,
            "Publication Date": pub_date_text,
            "Non-academic Author(s)": "; ".join(non_academic_authors),
            "Company Affiliation(s)": "; ".join(company_affiliations),
            "Corresponding Author Email": corresponding_email or "Not found"
        })

    return papers
