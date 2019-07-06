# market-color-map
This SW creates and display color-map of the historical market prices correlation:

Copy files to a local folder. Use your Web Browser (Google Chrome) to open SGraph.html.
Select csv file with historical Open, High, Low, Close data (as example you can use existing TEST.csv, or use Yahoo (https://help.yahoo.com/kb/SLN2311.html).
SGraph calculates correlation and potential profit (for one market share) by using correlation, historical data and the following strategy: Set "sell limit" every morning to "Price at Open" x "Calculated recommendation".

The X, Y color of colormap calculated on the basis of Time (X) and time difference between correlation intervals (Y).

These screen shots show Correlation colormap and graphic (caculate on the bases of the maximum correlation).

![Correlation colormap](/Test.png?raw=true)
