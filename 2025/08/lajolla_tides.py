import requests
from lxml import etree
from datetime import datetime

# NOAA tide prediction URL for La Jolla station (San Diego area)
TIDE_XML_URL = "https://tidesandcurrents.noaa.gov/noaatidepredictions.xml?station=9410230&datum=MLLW&units=english&time_zone=lst_ldt&interval=hilo&format=xml"

def fetch_low_tides(url):
    response = requests.get(url)
    response.raise_for_status()

    root = etree.fromstring(response.content)
    low_tides = []

    for item in root.xpath("//item"):
        event_type = item.xpath("eventType")[0].text
        if event_type.lower() == "low":
            date_time = item.xpath("date")[0].text
            height = item.xpath("height")[0].text
            low_tides.append((date_time, height))

    return low_tides

def display_tides(tides):
    print("📉 Upcoming Low Tides in San Diego (La Jolla Station):")
    for dt_str, height in tides:
        dt = datetime.strptime(dt_str, "%Y-%m-%d %H:%M")
        print(f"• {dt.strftime('%A, %b %d %I:%M %p')} — {height} ft")

if __name__ == "__main__":
    tides = fetch_low_tides(TIDE_XML_URL)
    display_tides(tides)

