import xml.etree.ElementTree as ET

my_xml = ET.parse('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/ShortDataset/JC_uncertain.xml' )
root = my_xml.getroot()

# print(root.tag)
#
# print(root.attrib)


# for id , child in enumerate(root):
#     # print(id)
#     print(child.attrib)
#     child.set('uncertain', "true" )


# print(ET.tostring(root, encoding='utf8').decode('utf8'))

# for seq in root.findall("./data"):
#     ET.SubElement(seq, 'sequence', )
#     seq.set('uncertain', "true" )
#     seq.attrib['value'] = "test"
#     print(seq.attrib['value'])
#
# ll = root.findall("./data")
# ET.SubElement(ll, 'sequence' )
data = root.find("data")
# el = ET.Element("XXXXXXXXXXXXXXXXXXXXXXXXXX")

for i in range(10):
    c = ET.Element("sequence")
    c.set("taxon" , str(i))
    c.set("uncertain" , "true")
    c.text = '\n' +  "ABCDEFG"+str(i) +'\n'
    data.insert(i,c)
    c.tail = "\n"


# new_dec = ET.SubElement(seq, "sequence")
# new_dec.set("taxon", "0")
# new_dec.set("uncertain", "true")

my_xml.write('/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/recombination/ShortDataset/test1.xml' ,encoding="utf-8", xml_declaration=True)

