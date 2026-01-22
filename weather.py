

import ee
from time import sleep
import csv
import datetime as dt
from newdatapts import sturgeon_pts, sacramento_pts, mississippi_pts,willamette_pts
service_account = 'justin@pajela.com'
credentials = ee.ServiceAccountCredentials(service_account, '/Users/justin/Downloads/hogue-428318-8f3882796b90.json')
ee.Initialize(credentials)


def getPtDataCloudCover(lon, lat, start, end):
    """Retrieve cloud cover data for a given point and time range."""
    print(f"Retrieving data for {lon}, {lat} from {start} to {end}")  
    region = ee.Geometry.Point(lon, lat)
    dataset = (
        ee.ImageCollection("NOAA/NCEP_DOE_RE2/total_cloud_coverage")
        .filterDate(start, end)
        .select('tcdc')  
        .mean()
        .clip(region)
    )
    
    mean_cloud_cover = dataset.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=5000,
        maxPixels=1e9
    )
    
    mean_cover_value = mean_cloud_cover.getInfo().get('tcdc')
    if mean_cover_value is not None:
        print(f"Cloud cover found: {mean_cover_value}")  
        return ee.Feature(region, {'mean_cloud_cover': mean_cover_value, 'date': start})
    else:
        print(f"No cloud cover data found for {start}")  
        return None

datapts = [
[[-88.4918473431673, 47.0322317215815], "Point 1"],
[[-88.4983930339682, 47.0291722945433], "Point 2"],
[[-88.500946237914, 47.023974996757], "Point 3"],
[[-88.5065058403125, 47.0172072368729], "Point 4"],
[[-88.5120662018321, 47.0093940825674], "Point 5"],
[[-88.5086818883236, 47.0018444871796], "Point 6"],
[[-88.5127251636897, 46.9934132623617], "Point 7"],
[[-88.5147371772906, 46.9887575049127], "Point 8"],
[[-88.5242113096453, 46.9841436922203], "Point 9"],
[[-88.5225062924311, 46.9774995349019], "Point 10"],
[[-88.5245877137213, 46.9703183109816], "Point 11"],
[[-88.5305066039486, 46.9630722786732], "Point 12"],
[[-88.5316716628541, 46.9581992667013], "Point 13"],
[[-88.5340780890439, 46.9497954834531], "Point 14"],
[[-88.5355263493911, 46.9458641655041], "Point 15"],
[[-88.5407419691544, 46.9374977627694], "Point 16"],
[[-88.5477193976257, 46.9291678162265], "Point 17"],
[[-88.5467595860259, 46.9220725451029], "Point 18"],
[[-88.538595512044, 46.9197851437687], "Point 19"],
[[-88.5301062334267, 46.9145114743785], "Point 20"],
[[-88.5204285516212, 46.9169583637674], "Point 21"],
[[-88.5221018399631, 46.9112035285609], "Point 22"],
[[-88.5213699573767, 46.9084199027455], "Point 23"],
[[-88.5217427991281, 46.9034476393227], "Point 24"],
[[-88.5231086281917, 46.9012489367373], "Point 25"],
[[-88.5257999545815, 46.8965227296531], "Point 26"],
[[-88.5286323043955, 46.889284778511], "Point 27"],
[[-88.5277859107066, 46.8869752783026], "Point 28"],
[[-88.5342569294217, 46.8817132222591], "Point 29"],
[[-88.5397764405922, 46.8779585289402], "Point 30"],
[[-88.5365850323907, 46.8720179000646], "Point 31"],
[[-88.5428580576994, 46.8647949711493], "Point 32"],
[[-88.5510349299235, 46.8663702576813], "Point 33"],
[[-88.5503454309566, 46.8598822873193], "Point 34"],
[[-88.5553305028551, 46.8532806251589], "Point 35"],
[[-88.5614529114476, 46.8528590726834], "Point 36"],
[[-88.5648233443727, 46.8498420678679], "Point 37"],
[[-88.571602762716, 46.8473670285299], "Point 38"],
[[-88.5751144065088, 46.8430831915474], "Point 39"],
[[-88.5818146780117, 46.8439266972668], "Point 40"],
[[-88.5881831024184, 46.8400704364522], "Point 41"],
[[-88.5966510444732, 46.834154922514], "Point 42"],
[[-88.6037727917441, 46.8281070008546], "Point 43"],
[[-88.6132361239808, 46.8261299633209], "Point 44"],
[[-88.6201705428804, 46.8263004921315], "Point 45"],
[[-88.624413511343, 46.819004170303], "Point 46"],
[[-88.6293915988674, 46.8158904312447], "Point 47"],
[[-88.6317783356078, 46.81144543956], "Point 48"],
[[-88.6291820318823, 46.804958221178], "Point 49"],
[[-88.6233724352426, 46.8003198398963], "Point 50"],
[[-88.621518871062, 46.7947292486894], "Point 51"],
[[-88.6229653584826, 46.7875411152177], "Point 52"],
[[-88.6177626474634, 46.783891990376], "Point 53"],
[[-88.6189109557813, 46.7777040373557], "Point 54"],
[[-88.618892105514, 46.7718099076622], "Point 55"],
[[-88.6190976774994, 46.7667636308083], "Point 56"],
[[-88.6173195100962, 46.7635320103913], "Point 57"],
[[-88.619328612475, 46.7561553062359], "Point 58"],
[[-88.6265588138926, 46.7547705485858], "Point 59"],
[[-88.6275795550997, 46.7517814603366], "Point 60"],
[[-88.6313918301411, 46.746426133719], "Point 61"],
[[-88.636045651259, 46.7417471973028], "Point 62"],
[[-88.64072096467, 46.7399129071823], "Point 63"],
[[-88.6421323069895, 46.7386947690867], "Point 64"],
[[-88.648107858947, 46.7347201239618], "Point 65"],
[[-88.6543976467458, 46.7329611134552], "Point 66"],
[[-88.6520863375024, 46.7288449746133], "Point 67"],
[[-88.6591068956028, 46.7280744292455], "Point 68"],
[[-88.6690070069148, 46.7257760653092], "Point 69"],
[[-88.6741669004285, 46.7179865276235], "Point 70"],
[[-88.6802972719124, 46.7098020894148], "Point 71"],
[[-88.6857857385092, 46.7011737523922], "Point 72"],
[[-88.6936833073728, 46.6939651889143], "Point 73"],
[[-88.7014926565683, 46.687924584729], "Point 74"],
[[-88.7112127752372, 46.6848335602869], "Point 75"],
[[-88.7155810870708, 46.6788067663126], "Point 76"],
[[-88.7203879205997, 46.672289374336], "Point 77"],
[[-88.7172066125891, 46.6679004802917], "Point 78"],
[[-88.7160058083944, 46.6628807383987], "Point 79"],
[[-88.7109264970739, 46.6607709435628], "Point 80"],
[[-88.7095886185288, 46.6529629163959], "Point 81"],
[[-88.7073277373229, 46.6496353130053], "Point 82"],
[[-88.6984874999338, 46.6470224101118], "Point 83"],
[[-88.6931538046895, 46.6423240887525], "Point 84"],
[[-88.6851041848428, 46.6383614159293], "Point 85"],
[[-88.6767533632346, 46.6370243169777], "Point 86"],
[[-88.6776969433631, 46.6294197693225], "Point 87"],
[[-88.679290173334, 46.6258890080235], "Point 88"],
[[-88.6789343504585, 46.6219333106775], "Point 89"],
[[-88.6760753977118, 46.6161454485521], "Point 90"],
[[-88.6736649149613, 46.6127038876564], "Point 91"],
[[-88.667817508961, 46.6052276714922], "Point 92"],
[[-88.662703405735, 46.5992933163632], "Point 93"],
[[-88.6698678097488, 46.5931296831827], "Point 94"],
[[-88.6666793082794, 46.5917431248101], "Point 95"],
[[-88.6746861153486, 46.5921702350387], "Point 96"],
[[-88.6763855525813, 46.5871881472605], "Point 97"],
[[-88.6781904905135, 46.5809555707538], "Point 98"],
[[-88.6693485276782, 46.5792206109084], "Point 99"],
[[-88.6695429115526, 46.5751663457497], "Point 100"]
]

times = [(str(year), str(day + 1),
          str(dt.date(year, 1, 1) + dt.timedelta(days=day)),
          str(dt.date(year, 1, 1) + dt.timedelta(days=day + 1)))
         for year in range(2020,2023)  # Adjust years as needed
         for day in range(365)]

def getCloudCoverTimeSeries(pts, times, output_file):
    all_results = []
    for pt in pts:
        lon, lat = pt[0]
        point_id = pt[1]
        print(f"Processing point: {point_id}, Longitude: {lon}, Latitude: {lat}")  
        for ts in times:
            year, day, start, end = ts
            try:
                feature = getPtDataCloudCover(lon, lat, start, end)
                if feature:
                    mean_cloud_cover = feature.getInfo()['properties']['mean_cloud_cover']
                    date = feature.getInfo()['properties']['date']
                    all_results.append({
                        'point_id': point_id, 
                        'lon': lon,  
                        'lat': lat,  
                        'date': date, 
                        'mean_cloud_cover': mean_cloud_cover
                    })
                    print(f"Data added for {point_id} on {date}") 
                else:
                    print(f"No data for {point_id} on {start}") 
            except Exception as e:
                print(f"Error retrieving data for {point_id} on {start}: {e}")

    with open(output_file, mode='w', newline='') as csvfile:
        fieldnames = ['point_id', 'lon', 'lat', 'date', 'mean_cloud_cover']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in all_results:
            writer.writerow(result)

if __name__ == "__main__":
    output_file = "cloud_cover_willametteriver.csv"
    getCloudCoverTimeSeries(willamette_pts, times, output_file)
    print(f"Data has been written to {output_file}")