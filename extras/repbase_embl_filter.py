# -*- coding: utf-8 -*-
import argparse

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-e", "--embl", help="EMBL file", required=True)

args = parser.parse_args()#pylint: disable=invalid-name

content_file = open(args.embl, 'r')
content = content_file.read()
content_split = content.split('//\n')
m_file = open('data/repbase_nonautonomous.csv','w') 
m_file.close()
for con in content_split:    
    file_id = con.strip().split('\n')[0].replace('ID','',1).strip().split(' ')[0]
    if file_id == '':
        continue
    if ('Nonautonomous' in con or 'non-autonomous' in con or 'Non-autonomous' in con) and \
        not 'Retrotransposon' in con:

        m_file = open('data/repbase_nonautonomous.csv','a') 
        m_file.write(file_id + '\n')
        m_file.close()