from obspy.io.xseed import Parser
import sys 

data1, data2 = sys.argv[1], sys.argv[2]

d1, d2 = Parser(data1), Parser(data2)

inv1, inv2 = d1.get_inventory(), d2.get_inventory()

for chn in inv1['channels']:
	channel_id, start_date, end_date, instrument = chn['channel_id'], chn['start_date'], chn['end_date'], chn['instrument']
	print channel_id
	PAZ = d1.get_paz(seed_id=channel_id,datetime=start_date)
	print PAZ

for chn in inv2['channels']:
	channel_id, start_date, end_date, instrument = chn['channel_id'], chn['start_date'], chn['end_date'], chn['instrument']
	print channel_id
	PAZ = d2.get_paz(seed_id=channel_id,datetime=start_date)
	print PAZ
