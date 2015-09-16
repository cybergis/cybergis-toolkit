rm test1.log
rm test1.out
rm test1.err
rm test1.nodelist

python test_nodes.py

sleep 60s

python parse_test1.py
