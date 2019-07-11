# market-color-map
This SW creates and display color-map of the historical market prices correlation:
Use WEB site https://andrei5055.github.io/ or 
Copy all files to a local folder and then open index.html with Google Chrome.
When you start it first time from local folder you need to use "Choose a file..." button and select local csv file with Open, High, Low and Close values (you can use provided TEST.csv file as example).
You can select your own csv file with historical Open, High, Low, Close data: use Yahoo (https://help.yahoo.com/kb/SLN2311.html) to download historical data.
SGraph calculates correlation and potential profit (for one market share) by using correlation, historical data and the following strategy: 
1. If you do not have shares, then buy when market open. 
2. Set "sell limit" to: "Price at Open" x "Calculated recommendation".

The X, Y color of colormap calculated on the basis of Time (X) and time difference between correlation intervals (Y). Red color - "good" correlation, Blue - bad.

These screen shots show Correlation colormap and profit graphic (caculate on the bases of the maximum correlation).

![Correlation colormap](/Test.png?raw=true)
