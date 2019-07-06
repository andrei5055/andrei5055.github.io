// SGraph GUI support
/*
Copyright (c) 2019, Leonid M Tertitski (Leonid.M.Tertitski@gmail.com)
All rights reserved.

Redistribution and use in source and binary forms, without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
function sgraphgui() {
    var MATRIX_Y = 200;
    var MAX_ROWS_IN_RESULT_TABLE = 100;

    var m_historicalDataParamsInd = 3;
    var m_nParamsWindows = 5; // number of sg_input_paramsX_div parameters windows in HTML
    var m_zIndexParamsStart = 11;
    var m_zIndexParamsMax = m_zIndexParamsStart + m_nParamsWindows - 1;

    var m_params_windows = [ // must be the m_nParamsWindows+2 values 
        "Show all",
        "Checkboxes",
        "Sliders",
        "Historical Data",
        "Archived Parameters",
        "Order for 'Tomorrow'",
        "Hide all"];
 
    var m_showParamWindow = [false, false, false, false, true]; // by default show Order for 'Tomorrow' only. Must be m_nParamsWindows values 

    var m_nSliders = 6;
    var xOffset = 0;
    var angleOffset = 0;
    var m_CorrLength = 16;
    var m_OrderScale = 0.1;
    var m_AverageLengthMax = 100;
    var m_MaxLost = 0.4;
    var m_AverageLength = 25;
    var m_Spread = 0.03;
    var m_bShowOpen = true;
    var m_bShowCandleStick = true;
    var m_bKnownOrder = true;
    var m_bHighlightCorrIntervals = true;
    var m_bShowAverage = false;
    var m_bShowVolume = false;
    var m_bShowPerfomance = false;
    var m_bShowOrder = false;
    var m_bShowProfit = false;
    var m_bShowCorr = false;
    var m_bShowCorrPos = false;
    var m_bUseLeoDraw = true;
    var m_bUseCanvas = true;
    var m_bCalculatedOrder = true;
    var m_bFixedOrder = false;
    var m_nDaysToKeepOpen = 1;
    var m_movingId = 0;
    var m_zoomSliderId = 0;
    var m_faviconIndex = 0;
    var m_paramsSliderTimer = 0;
    var m_firstTime = true;
    var m_numberOfColumnsInResult = -1;
    var m_profit = 0.0;
    var m_archiveTableChangedIndex = -1;
    var m_mouseXpos = 0;
    var m_mouseYpos = 0;
    var m_cnv = null;
    var m_ctx = null;
    var m_sgraph = null;
    var m_lastClickTime = new Date();
    var m_lastArchiveClickTime = new Date();
    var m_lastPasteTime = new Date();
    var m_TouchScreen = false;
    var m_mouseupEvent = 'mouseup';
    var m_mousemoveEvent = 'mousemove';
    var m_mousedownEvent = 'mousedown';
    var m_copyReadyTm = 0;
    function onLoadSetup() {
        let w = window,
            d = document,
            e = d.documentElement,
            g = d.getElementsByTagName('body')[0],
            x = w.innerWidth || e.clientWidth || g.clientWidth,
            y = w.innerHeight || e.clientHeight || g.clientHeight;
        m_TouchScreen = 'ontouchstart' in window;
        m_mouseupEvent = m_TouchScreen ? 'touchend' : 'mouseup';
        m_mousemoveEvent = m_TouchScreen ? 'touchmove' : 'mousemove';
        m_mousedownEvent = m_TouchScreen ? 'touchstart' : 'mousedown';
        m_cnv = document.getElementById('sg_canvas');
        m_ctx = m_cnv.getContext('2d');

        // add mouse click listener to graphic area
        if (m_cnv) {
            m_cnv.addEventListener("click", mouseOnGraphClick, false);
        }
        // Setup params windows names
        setupParamsWindows();

        // init Start/Stop for 'zoom slider'
        zoomSliderInit();

        // set mouse click calbacks for parameters windows
        setParamsWindowsClick();

        getSavedParams();

        // init m_sgraph
        m_sgraph = sgraph();
        m_sgraph.onDrawEnd = sgraphCalculationEnd;
        m_sgraph.onOptimizeEnd = onOptimizationEnd;
        m_sgraph.setup(x, y);

        // hide start/stp zoom slider
        hideShowSlider();
        document.addEventListener(m_mousedownEvent, processWindowOnClick, false);
        // setup dynamic favicon
        drawFavicon();

        if (localStorage && m_sgraph) {
            let strHistData = localStorage.getItem("sg_hist_data");
            let strHistName = localStorage.getItem("sg_hist_name");
            if (!strHistData || !strHistName) {
                loadFile("TEST.csv");              
            }
            else {
                document.getElementById("sg_filename").innerText = strHistName;
                m_sgraph.histDataToTable(strHistData);
                showSGraph();
            }
        }
    }
    function loadFile(filePath) {
        var xmlhttp = new XMLHttpRequest();
        xmlhttp.onload = historyDataLoaded;
        xmlhttp.open("GET", filePath, true);
        xmlhttp.setRequestHeader("Content-Type", "text/html");
        xmlhttp.send();
    }
    function historyDataLoaded() {
        if (this.status == 200 && this.responseText) {
            document.getElementById("sg_filename").innerText = "TEST.csv";
            m_sgraph.histDataToTable(this.responseText);
            showSGraph();
        } else {
            document.body.innerHTML = this.responseText;
        }
    }

    function getMouseX(e) {
        return (m_TouchScreen && e.touches) ? e.touches[0].clientX : e.clientX;
    }
    function getMouseY(e) {
        return (m_TouchScreen && e.touches) ? e.touches[0].clientY : e.clientY;
    }
    function paramsWindowClick(e) {
        e = e || window.event;
        //e.stopPropagation();
        pid = e.target.id;
        if (pid === "")
            return;
        let tm = new Date();
        let timeSinceLastClick = tm.getTime() - m_lastClickTime.getTime();
        //console.log(timeSinceLastClick);
        if (timeSinceLastClick < 300) {
            let dt = m_lastClickTime.getTime() - m_lastArchiveClickTime.getTime();
            //console.log("dt="+dt);
            if (dt > 500)
                hideParamsWindow(pid);
        }
        else {
            m_lastClickTime = new Date();
            paramsWindowToTop(pid);
        }
    }
    function paramsWindowToTop(pid) {
        if (isNaN(pid)) {
            if (!pid.startsWith("sg_input"))
                return;
            ind = pid.replace(/\D/g, "");
        }
        else
            ind = pid;
        let paramsTop = getParamsWindowByInd(ind);
        if (paramsTop) {
            paramsTop.style.display = "inline-block";
            let zind = parseInt(paramsTop.style.zIndex);
            m_showParamWindow[ind - 1] = true;
            if (zind != m_zIndexParamsMax) {
                for (let i = 1; i <= m_nParamsWindows; i++) {
                    let params = getParamsWindowByInd(i);
                    if (params && parseInt(params.style.zIndex) > zind) {
                        params.style.zIndex = parseInt(params.style.zIndex) - 1;
                        saveParamsWindowPos(ind);
                    }
                }
                paramsTop.style.zIndex = m_zIndexParamsMax;
            }
            saveParamsWindowPos(ind);
        }
    }
    function dragMouseDown(e) {
        e = e || window.event;
        e.preventDefault();
        // get the mouse cursor position at startup:
        m_mouseXpos = getMouseX(e);
        m_mouseYpos = getMouseY(e);

        m_movingId = e.target.id;
        paramsWindowToTop(m_movingId);

        document.addEventListener(m_mouseupEvent, closeDragElement, false);
        document.addEventListener(m_mousemoveEvent, elementDrag, false);
    }
    function setParamsWindowPos(params, newX, newY) {
        let minwidth = 20; // params.offsetWidth;
        let minHeight = 20; // params.offsetHeight;
        if (newX + params.offsetWidth < minwidth)
            newX = minwidth - params.offsetWidth;
        else if (newX + minwidth > window.innerWidth)
            newX = window.innerWidth - minwidth;
        params.style.left = newX + "px";

        if (newY < 0)
            newY = 0;
        else if (newY + minHeight > window.innerHeight)
            newY = window.innerHeight - minHeight;
        params.style.top = newY + "px";
    }
    function elementDrag(e) {
        e = e || window.event;
        //e.preventDefault();
        // calculate the new cursor position:
        let posx = m_mouseXpos - getMouseX(e);
        let posy = m_mouseYpos - getMouseY(e);
        m_mouseXpos = getMouseX(e);
        m_mouseYpos = getMouseY(e);
        // set the element's new position:
        let params = document.getElementById(m_movingId + "_div");
        if (params) {
            let newX = params.offsetLeft - posx;
            let newY = params.offsetTop - posy;
            setParamsWindowPos(params, newX, newY);
        }
    }

    function closeDragElement() {
        document.removeEventListener(m_mousedownEvent, closeDragElement, false);
        document.removeEventListener(m_mousemoveEvent, elementDrag, false);
        saveParamsWindowsPos();
    }

    function zoomSliderInit() {
        let slider = document.getElementById("zoom_slider_div");
        let bar = document.getElementById("zoom_bar");
        let start = document.getElementById("zoom_bar_start");
        let stop = document.getElementById("zoom_bar_stop");
        if (m_TouchScreen) {

            slider.ontouchstart = zoomSliderMouseDown;
            bar.ontouchstart = zoomSliderMouseDown;
            start.ontouchstart = zoomSliderMouseDown;
            stop.ontouchstart = zoomSliderMouseDown;
        }
        else {
            slider.onmousedown = zoomSliderMouseDown;
            bar.onmousedown = zoomSliderMouseDown;
            start.onmousedown = zoomSliderMouseDown;
            stop.onmousedown = zoomSliderMouseDown;
        }
    }
    function zoomSliderMouseDown(e) {
        e = e || window.event;
        e.preventDefault();
        // get the mouse cursor position at startup:
        m_mouseXpos = getMouseX(e);
        m_zoomSliderId = e.target.id;
        document.addEventListener(m_mouseupEvent, closezoomSlider, false);
        document.addEventListener(m_mousemoveEvent, zoomSliderDrag, false);
    }
    function zoomSliderDrag(e) {
        e = e || window.event;
        //e.preventDefault();
        // calculate the new cursor position:
        let pos1 = m_mouseXpos - getMouseX(e);
        m_mouseXpos = getMouseX(e);
        // set the element's new position:
        let bar = document.getElementById("zoom_bar");
        let slider = document.getElementById("zoom_slider_div");
        let barStart = document.getElementById("zoom_bar_start");
        let barStop = document.getElementById("zoom_bar_stop");
        if (bar && slider && barStop && barStop && m_sgraph && m_sgraph.m_nData > 1 && m_zoomSliderId) {
            let sliderWidth = slider.offsetWidth;
            let minScrolBarWidth = barStart.offsetWidth + barStop.offsetWidth;
            if (sliderWidth >= minScrolBarWidth) {
                let x0 = m_sgraph.m_Start * sliderWidth / m_sgraph.m_nData;
                let x1 = m_sgraph.m_Stop * sliderWidth / m_sgraph.m_nData;
                let dx = x1 - x0;
                if (m_zoomSliderId == "zoom_bar_start") {
                    x0 = x0 - pos1;
                    if (x0 < 0)
                        x0 = 0;
                    else if (x0 > x1 - minScrolBarWidth)
                        x0 = x1 - minScrolBarWidth;
                }
                else if (m_zoomSliderId == "zoom_bar_stop") {
                    x1 = x1 - pos1;
                    if (x1 < x0 + minScrolBarWidth)
                        x1 = x0 + minScrolBarWidth;
                    else if (x1 > sliderWidth)
                        x1 = sliderWidth;
                }
                else {
                    x0 = x0 - pos1;
                    x1 = x1 - pos1;
                    if (x0 < 0) {
                        x0 = 0;
                        x1 = dx;
                    }
                    else if (x1 > sliderWidth) {
                        x1 = sliderWidth;
                        x0 = x1 - dx;
                    }
                }
                setGraphStartStop(x0 / sliderWidth, x1 / sliderWidth);
            }
        }
    }
    function zoomSliderSetStartStop() {

        let bar = document.getElementById("zoom_bar");
        let start = document.getElementById("zoom_bar_start");
        let stop = document.getElementById("zoom_bar_stop");
        if (bar && m_sgraph && start && stop) {
            let rect = bar.getBoundingClientRect();
            start.style.left = rect.left + "px";
            stop.style.left = (rect.right - stop.offsetWidth) + "px";
        }
    }

    function closezoomSlider() {
        document.removeEventListener(m_mouseupEvent, closezoomSlider, false);
        document.removeEventListener(m_mousemoveEvent, zoomSliderDrag, false);
        saveParamsValues();
    }
    function setParameters() {
        m_bShowOpen = document.getElementById("param_cb21").checked;
        m_bShowCandleStick = document.getElementById("param_cb22").checked;
        m_bKnownOrder = document.getElementById("param_cb23").checked;
        m_bHighlightCorrIntervals = document.getElementById("param_cb24").checked;
        m_bShowAverage = document.getElementById("param_cb25").checked;
        m_bShowVolume = document.getElementById("param_cb26").checked;
        m_bShowOrder = document.getElementById("param_cb27").checked;
        m_bShowProfit = document.getElementById("param_cb28").checked;
        m_bShowCorr = document.getElementById("param_cb29").checked;
        m_bShowCorrPos = document.getElementById("param_cb30").checked;
        m_bFixedOrder = document.getElementById("param_cb31").checked;
        m_bCalculatedOrder = document.getElementById("param_cb33").checked;
        m_CorrLength = parseInt(document.getElementById("param_sl1").value);
        m_OrderScale = parseFloat(document.getElementById("param_sl2").value);
        m_AverageLength = parseInt(document.getElementById("param_sl3").value);
        m_AverageLengthMax = parseFloat(document.getElementById("param_sl3").max);
        m_MaxLost = parseFloat(document.getElementById("param_sl4").value);
        m_nDaysToKeepOpen = parseInt(document.getElementById("param_sl5").value);
        m_Spread = parseFloat(document.getElementById("param_sl6").value);
        updateSlidersLabels();
        if (m_bKnownOrder && m_bFixedOrder) {
            m_bFixedOrder = false;
            document.getElementById("param_cb31").checked = false;
        }
        showVisibleParamsWindows();
    }
    function resizeGraph() {
        let w = window,
                d = document,
                e = d.documentElement,
                g = d.getElementsByTagName('body')[0],
                x = w.innerWidth || e.clientWidth || g.clientWidth,
                y = w.innerHeight || e.clientHeight || g.clientHeight;
        m_cnv = document.getElementById('sg_canvas');
        m_ctx = m_cnv.getContext('2d');
        m_cnv.width = x;
        m_cnv.height = y - m_cnv.offsetTop;
        if (m_sgraph) {
            m_sgraph.setup(x, y);
        }
        showSGraph();
    }
    function redrawGraph() {
        showSGraph();
    }
    function onZoomIn() {
        if (m_sgraph) {
            m_sgraph.zoom(1);
            m_sgraph.processLoadedData();
        }
    }
    function onZoomOut() {
        if (m_sgraph) {
            m_sgraph.zoom(-1);
            m_sgraph.processLoadedData();
        }
    }
    function tableHeaderItem(columnName) {
        return "<th>" + columnName + "</th>";
    }
    function tableCellItem(cellValue) {
        return "<td>" + cellValue + "</td>";
    }
    function createArhiveTableHead() {
        let table = "<table id='sg_params_table'>";
        table += "<thead class='theadTable'>";
        table += "<tr>";
        table += tableHeaderItem("Label");
        table += tableHeaderItem("File Name");
        table += tableHeaderItem("Start");
        table += tableHeaderItem("End");
        m_nSliders = 0;
        for (m_nSliders = 0; m_nSliders < 20; m_nSliders++) {
            let vname = document.getElementById("sg_param_label" + (m_nSliders + 1));
            if (!vname)
                break;
            let vnameh = vname.textContent.split(":");
            table += tableHeaderItem(vnameh[0]);
        }
        table += tableHeaderItem("Mode");
        table += tableHeaderItem("Profit($)");
        table += "</tr>";
        table += "</thead>";
        table += "<tbody>";
        return table;
    }
    function createArchiveTableTail(table) {
        table += '</tbody>';
        table += '</table>';
        return table;
    }
    function getSavedParams() {
        let paramsTable = document.getElementById("sg_params_tbl");
        let zTotal = 0;
        for (let i = 1; i <= m_nParamsWindows; i++) {
            let params = getParamsWindowByInd(i);
            if (params) {
                params.addEventListener(m_mousedownEvent, sg_main.paramsWindowClick, false);
                if (localStorage) {
                    let x = localStorage.getItem("sg_params_window" + i + "_x");
                    let y = localStorage.getItem("sg_params_window" + i + "_y");
                    let z = parseInt(localStorage.getItem("sg_params_window" + i + "_z"));
                    let show = localStorage.getItem("sg_params_window" + i + "_show");
                    if (show)
                        m_showParamWindow[i - 1] = (show === "show") ? true : false;
                    let xNum = parseFloat(x, 10);
                    let yNum = parseFloat(y, 10);
                    if (!isNaN(xNum) && !isNaN(yNum)) {
                        setParamsWindowPos(params, xNum, yNum);
                    }
                    if (z) {
                        params.style.zIndex = z;
                        zTotal += z;
                    }
                }
            }
        }
        if (zTotal != (m_zIndexParamsMax + m_zIndexParamsStart) * (m_zIndexParamsMax - m_zIndexParamsStart + 1) / 2) {
            for (let i = 1; i <= m_nParamsWindows; i++) {
                let params = getParamsWindowByInd(i);
                if (params)
                    params.style.zIndex = m_zIndexParamsStart + i - 1;
            }
        }
        if (localStorage) {
            let ids = document.querySelectorAll('[id]');

            Array.prototype.forEach.call(ids, function (el, i) {
                if (el.id.startsWith('param')) {
                    let vsaved = localStorage.getItem(el.id);
                    if (vsaved !== null) {
                        let val;
                        if (el.type == "checkbox" || el.type == "radio")
                            el.checked = vsaved === "true";
                        else if (el.type != "file")
                            el.value = vsaved;
                    }
                }
            });
        }
        if (paramsTable) {
            let table = createArhiveTableHead();
            if (m_nSliders > 0) {
                m_numberOfColumnsInResult = m_nSliders + 6;
                if (localStorage) {
                    let touchEvent = m_TouchScreen ? 'ontouchstart' : 'onmousedown';
                    for (let i = 1; i <= MAX_ROWS_IN_RESULT_TABLE; i++) {
                        let rowTxt = localStorage.getItem("sg_params_table" + i);
                        if (!rowTxt || rowTxt.startsWith("undefined") || rowTxt.startsWith("0,1,"))
                            break;
                        let rowCells = rowTxt.split(",");
                        table += "<tr " + touchEvent + "='sg_main.resutlTableRowClicked(this)'>";
                        if (m_numberOfColumnsInResult == rowCells.length) {
                            for (let rowCell = 0; rowCell < rowCells.length; rowCell++) {
                                if (rowCell == 0)
                                    table += "<td tabindex='" + i + "' contenteditable='true' spellcheck = 'false' style='user-select: text' onblur='sg_main.storeTableIfNeeded()' oninput='sg_main.tableChanged(this)'>" + rowCells[rowCell] + "</td>";
                                else
                                    table += tableCellItem(rowCells[rowCell]);
                            }
                        }
                        table += "</tr>";
                    }
                }
            }
            table = createArchiveTableTail(table);
            paramsTable.innerHTML = table;
        }
        setParameters();
    }
    function setResultString() {
        let resultDiv = document.getElementById("sg_order_tomorrow");
        if (m_sgraph && resultDiv) {
            resultDiv.style = "background-color: #b0ffff";
            resultDiv.innerHTML = m_sgraph.orderString("sg_order_tomorrow");
        }

        let result = document.getElementById("sg_result");
        let fileName = document.getElementById("sg_filename");
        if (result && fileName && m_sgraph && m_sgraph.m_nData > 0) {
            if (localStorage != undefined) {
                localStorage.setItem("sg_hist_name", fileName.innerHTML);
            }
            let resultTxt = fileName.innerText;
            let iStart = Math.floor(m_sgraph.m_Start);
            let iStop = Math.floor(m_sgraph.m_Stop);
            if (iStart < 0)
                iStart = 0;
            if (iStop >= m_sgraph.m_nData)
                iStop = m_sgraph.m_nData - 1;
            resultTxt += "," + m_sgraph.m_DateTime[iStart] + "," + m_sgraph.m_DateTime[iStop];
            for (let i = 1; i <= m_numberOfColumnsInResult; i++) {
                let vslider = document.getElementById("param_sl" + i);
                if (!vslider)
                    break;
                resultTxt += "," + vslider.value;
            }
            let mode = m_bFixedOrder ? "Fixed order" : m_bKnownOrder ? "Known price" : "Estimated order";
            resultTxt += "," + mode;
            resultTxt += "," + m_profit.toString();
            result.innerText = resultTxt;
        }
    }
    function copyParamsUsedToTable(strToCopy) {
        let table = document.getElementById("sg_params_table");
        if (table && strToCopy) {
            let strCopy = strToCopy + '\n';
            let allRows = strToCopy.split(/\r?\n|\r/);
            let rowIndex = 1;
            for (let singleRow = 0; singleRow < allRows.length; singleRow++) {
                let vrow = allRows[singleRow];
                if (vrow.startsWith("Label,File Name,") || vrow.startsWith("Label\tFile Name\t"))
                    continue; // ignore header
                let vcells = vrow.split(/,|\t/);
                if (vcells.length != m_numberOfColumnsInResult) {
                    if (vrow.length > 0)
                        alert("Error in line " + (singleRow + 1) + ":\n     " + vrow + "\n\n" + "Number of comma separated values " + vcells.length + ", expected " + m_numberOfColumnsInResult);
                    break;
                }
                let row = table.insertRow(rowIndex);
                row.addEventListener(m_mousedownEvent, function () { resutlTableRowClicked(this); }, false);
                rowIndex++;
                for (let i = 0; i < vcells.length; i++) {
                    let cell = row.insertCell(i);
                    if (i == 0) {
                        cell.tabindex = singleRow+1;
                        cell.contentEditable = true;
                        cell.style.userSelect = 'text';
                        cell.spellcheck = false;
                        //cell.addEventListener("blur", storeTableIfNeeded, false);
                        //cell.addEventListener("input", function () { tableChanged(this); }, false);
                        cell.onblur = storeTableIfNeeded;
                        cell.oninput = function () { tableChanged(this); };
                    }
                    cell.innerText = vcells[i];
                }
                for (let i = MAX_ROWS_IN_RESULT_TABLE + 1; i < table.rows.length; i++) {
                    table.deleteRow(i);
                }
            }
            saveParamsTable(-1);
        }
    }
    function saveParamsWindowPos(ind) {
        let z = 0;
        if (localStorage != undefined) {
            let params = getParamsWindowByInd(ind);
            if (params) {
                let rect = params.getBoundingClientRect();
                z = parseInt(params.style.zIndex);
                if (rect && rect.width != 0 && rect.height != 0) {
                    setParamsWindowPos(params, params.offsetLeft, params.offsetTop);
                    localStorage.setItem("sg_params_window" + ind + "_x", params.offsetLeft.toString());
                    localStorage.setItem("sg_params_window" + ind + "_y", params.offsetTop.toString());
                }
                localStorage.setItem("sg_params_window" + ind + "_show", (m_showParamWindow[ind - 1] === true) ? "show" : "hide");
            }
        }
        return z;
    }
    function saveParamsWindowsPos() {
        if (localStorage != undefined) {
            for (let i = 1; i <= m_nParamsWindows; i++) {
                saveParamsWindowPos(i);
            }
        }
    }
    function saveParamsValues() {
        if (localStorage != undefined) {
            let ids = document.querySelectorAll('[id]');
            Array.prototype.forEach.call(ids, function (el, i) {
                if (el.id.startsWith('param')) {
                    let val;
                    if (el.type == "checkbox" || el.type == "radio")
                        val = el.checked;
                    else
                        val = el.value;
                    localStorage.setItem(el.id, val);
                }
            });
        }
    }
    function saveParamsTable(index) {
        if (localStorage != undefined) {
            let table = document.getElementById("sg_params_table");
            if (table && table.rows && table.rows.length > 0) {
                let trows = table.rows;
                let ind1 = (index < 1) ? 1 : index;
                let ind2 = (index < 1) ? trows.length - 1 : index;
                for (let i = ind1; i <= ind2; i++) {
                    let vrow = getTableRowStr(table, i);
                    localStorage.setItem("sg_params_table" + i, vrow);
                }
            }
        }
    }
    function saveParamsAll() { // Do we need it?
        saveParamsWindowsPos();
        saveParamsValues();
        saveParamsTable(-1);
    }
    function copyTableToString() {
        let table = document.getElementById("sg_params_table");
        let tableStr = "";
        if (table && table.rows && table.rows.length > 0) {
            let trows = table.rows;
            for (let i = 0; i < trows.length; i++) {
                tableStr += getTableRowStr(table, i) + '\n';
            }
        }
        return tableStr;
    }
    function copyTableToClipboard() {
        let str = copyTableToString();
        copyStringToClipboard(str);
    }
    function getTableRowStr(table, iRow) {
        let row = table.rows[iRow];
        let oCells = row.cells;
        let vrow = "";
        // gets amount of cells of current row
        let cellLength = oCells.length;
        // loops through each cell in current row
        for (let j = 0; j < cellLength; j++) {
            if (j > 0)
                vrow += ",";
            vrow += oCells[j].innerText;
        }
        return vrow;
    }
    function showSGraph() {
        setParameters();
        saveParamsValues();
        document.querySelector("body").style['overflow-y'] = 'hidden';
        if (m_sgraph) {
            m_sgraph.processLoadedData();
        }
    }
    function showNewSGraph() {
        if (m_sgraph) {
            showSGraph();
        }
    }
    function clickArhiveTable() {
        let copyAttr = document.getElementById("sg_paste_table");
        if (copyAttr) {
            copyAttr.select();
        }
    }
    function copyPasteArchiveFinished() {
        clearTimeout(m_copyReadyTm);
        let copyAttr = document.getElementById("sg_paste_table");
        if (copyAttr) {
            copyAttr.value = "Copy/Paste Table";
            copyAttr.select();
        }
    }
    function copyArchiveReady() {
        clearTimeout(m_copyReadyTm);
        let copyAttr = document.getElementById("sg_paste_table");
        if (copyAttr) {
            copyAttr.value = "Done";
            m_copyReadyTm = window.setTimeout(copyPasteArchiveFinished, 1000);
        }
    }
    function copyArhiveTable() {
        let copyAttr = document.getElementById("sg_paste_table");
        if (copyAttr) {
            copyAttr.value = copyTableToString();
            copyAttr.select();
            m_copyReadyTm = window.setTimeout(copyArchiveReady, 50);
        }
    }
    function pasteArhiveTable() {
        let pasteAttr = document.getElementById("sg_paste_table");
        if (pasteAttr) {
            pasteAttr.value = "";
            m_lastPasteTime = new Date();
        }
    }
    function rejectArhiveInput() {
        let pasteAttr = document.getElementById("sg_paste_table");
        if (pasteAttr) {
            let tm = new Date();
            if (tm.getTime() - m_lastPasteTime.getTime() < 300) {
                copyParamsUsedToTable(pasteAttr.value);
                pasteAttr.value = "Done";
                m_copyReadyTm = window.setTimeout(copyPasteArchiveFinished, 1000);
            }
            else
                m_copyReadyTm = window.setTimeout(copyPasteArchiveFinished, 10);
        }
    }
    function copyResultToTable() {
        let resultAttr = document.getElementById("sg_result");
        if (resultAttr) {
            copyParamsUsedToTable("," + resultAttr.innerHTML);
        }
    }
    function copyStringToClipboard(str) {
        // Create new element
        let input = document.createElement('textarea');
        input.setAttribute("readonly", "");
        input.style = { position: 'absolute', left: '-9999px' };
        input.value = str;
        document.body.appendChild(input);
        // Select text inside element
        input.select();
        // Copy text to clipboard
        document.execCommand("copy");
        document.body.removeChild(input);
    }
    function setParamsWindowsClick() {
        for (let i = 1; i <= m_nParamsWindows; i++) {
            let paramsHeader = document.getElementById("sg_input_params" + i);
            if (paramsHeader) {
                if (m_TouchScreen)
                    paramsHeader.ontouchstart = dragMouseDown;
                else
                    paramsHeader.onmousedown = dragMouseDown;
            }
        }
    }
    function setupParamsWindows() {
        let pd_menu = document.getElementsByName("sg_pd_menu");
        let ind = 0;
        for (let i = 0; i < pd_menu.length; i++, ind++) {
            if (ind < m_params_windows.length) {
                pd_menu[i].innerText = m_params_windows[ind];
            }
        }
        for (let i = 1; i <= m_nParamsWindows && i < m_params_windows.length; i++) {
            let paramsHeader = document.getElementById("sg_input_params" + i);
            if (paramsHeader)
                paramsHeader.innerText = m_params_windows[i];
        }
    }
    function handleFile(files) {
        if (m_sgraph && files && files.length > 0) {
            document.getElementById("sg_filename").innerText = files[0].name;
            m_sgraph.loadData(files[0]);
        }
    }
    function showVisibleParamsWindows() {
        for (let i = 1; i <= m_nParamsWindows; i++) {
            let params = getParamsWindowByInd(i);
            params.style.display = (m_showParamWindow[i - 1] && m_sgraph && m_sgraph.m_nData > 0) ? "inline-block" : "none";
        }
    }
    function showHideAllParamsWindows(show) {
        for (let i = 1; i <= m_nParamsWindows; i++) {
            let params = getParamsWindowByInd(i);
            params.style.display = (show && m_sgraph && m_sgraph.m_nData > 0) ? "inline-block" : "none";
            m_showParamWindow[i - 1] = params.style.display != "none";
        }
        saveParamsWindowsPos();
    }
    function drawFavicon() {
        let cnv = document.createElement('canvas');
        let fsize = 192;
        cnv.width = fsize; cnv.height = fsize;
        let lw = 3;
        let iw = fsize;
        let ih = fsize;
        let ix = (fsize - iw) / 2 + lw / 2;
        let iy = lw / 2;
        let corner = fsize / 6;
        let ctx = cnv.getContext('2d');
        ctx.imageSmoothingEnabled = true;
        ctx.fillStyle = "#FFF";
        ctx.fillRect(0, 0, cnv.width, cnv.height);
        ctx.strokeStyle = "#eee";
        ctx.fillStyle = "#eee";
        ctx.lineWidth = lw;
        // draw rectangle with "round" corners
        ctx.beginPath();
        ctx.lineTo(ix + corner, iy);
        ctx.lineTo(fsize - ix - corner, iy);
        ctx.arc(fsize - ix - corner, iy + corner, corner, 1.5 * Math.PI, 2 * Math.PI);
        ctx.lineTo(fsize - ix, fsize - iy - corner);
        ctx.arc(fsize - ix - corner, fsize - iy - corner, corner, 0, 0.5 * Math.PI);
        ctx.lineTo(ix + corner, fsize - iy);
        ctx.arc(ix + corner, fsize - iy - corner, corner, 0.5 * Math.PI, Math.PI);
        ctx.lineTo(ix, iy + corner);
        ctx.arc(ix + corner, iy + corner, corner, Math.PI, 1.5 * Math.PI);
        ctx.closePath();
        ctx.stroke();
        ctx.fill();
        ctx.strokeStyle = "#00f";
        ctx.fillStyle = "#00f";
        let r = fsize / 16;
        switch (m_faviconIndex) {
            case 0:
                // Draw "3"
                for (let i = 0; i < 3; i++) {
                    ctx.beginPath();
                    ctx.arc(fsize / 2, fsize * (i + 1) / 4, r, 0, 2 * Math.PI);
                    ctx.fill();
                }
                break;
            case 1:
                // Draw "7"
                for (let i = 0; i < 3; i++) {
                    ctx.beginPath();
                    ctx.arc(fsize / 4, fsize * (i + 1) / 4, r, 0, 2 * Math.PI);
                    ctx.fill();
                    ctx.beginPath();
                    ctx.arc(fsize * 3 / 4, fsize * (i + 1) / 4, r, 0, 2 * Math.PI);
                    ctx.fill();
                }
                ctx.beginPath();
                ctx.arc(fsize / 2, fsize * 3 / 8, r, 0, 2 * Math.PI);
                ctx.fill();
                break;
            case 2:
                // Draw ace
                ctx.beginPath();
                ctx.arc(ix + lw + r * 2, lw + r * 2, r, 0, 2 * Math.PI);
                ctx.fill();
                ctx.beginPath();
                ctx.arc(iw - lw - r * 2, ih - lw - r * 2, r, 0, 2 * Math.PI);
                ctx.fill();
                let ra = r * 2;
                ctx.beginPath();
                ctx.arc(fsize / 2, fsize / 2 - ra - r, ra, 0, 2 * Math.PI);
                ctx.fill();
                ctx.beginPath();
                ctx.arc(fsize / 2 - ra, fsize / 2 - r, ra, 0, 2 * Math.PI);
                ctx.fill();
                ctx.beginPath();
                ctx.arc(fsize / 2 + ra, fsize / 2 - r, ra, 0, 2 * Math.PI);
                ctx.fill();
                ctx.lineWidth = 1;
                ctx.beginPath();
                ctx.moveTo(fsize / 2, fsize / 2);
                ctx.lineTo(fsize / 2 - ra, fsize / 2 + ra * 2);
                ctx.lineTo(fsize / 2 + ra, fsize / 2 + ra * 2);
                ctx.closePath();
                ctx.stroke();
                ctx.fill();
                break;
        }
        m_faviconIndex = (m_faviconIndex + 1) % 3;
        let link = document.createElement('link');
        link.type = 'image/x-icon';
        link.rel = 'shortcut icon';
        link.href = cnv.toDataURL("image/x-icon");
        let h = document.getElementsByTagName('head')[0];
        if (h.hasChildNodes())
            h.replaceChild(link, h.childNodes[0]);
        else
            h.appendChild(link);
        window.setTimeout(drawFavicon, 2000);
    }
    function sgraphCalculationEnd() {
        setSliderStartStop();
        setResultString(); // must be after setSliderStartStop
        showVisibleParamsWindows();
        m_copyReadyTm = window.setTimeout(copyPasteArchiveFinished, 100);
    }
    function hideShowSlider() {
        let dateSlider = document.getElementById("zoom_slider");
        if (!m_sgraph || m_sgraph.m_nData < 1)
            dateSlider.style.display = "none";
        else
            dateSlider.style.display = "inline-block";
        return true;
    }
    function setSliderStartStop() {
        let canvasDiv = document.getElementById("div_canvas");
        let canvas = document.getElementById("sg_canvas");
        let sliderDiv = document.getElementById("zoom_slider");
        let dateStart = document.getElementById("zoom_slider_start_date");
        let dateStop = document.getElementById("zoom_slider_stop_date");
        let slider = document.getElementById("zoom_slider_div");
        let bar = document.getElementById("zoom_bar");
        let barStart = document.getElementById("zoom_bar_start");
        let barStop = document.getElementById("zoom_bar_stop");
        let historicalData = document.getElementById("sg_input_params" + m_historicalDataParamsInd);
        let fileName = document.getElementById("sg_filename");
        if (hideShowSlider() && dateStart && dateStop && slider && bar && barStart && barStop && historicalData && fileName) {
            historicalData.textContent = m_params_windows[m_historicalDataParamsInd] + ": " + fileName.innerText;
            let sliderWidth = canvas.offsetWidth - dateStart.offsetWidth - dateStop.offsetWidth;
            sliderDiv.style.left = canvas.offsetLeft + "px";
            sliderDiv.style.width = canvas.offsetWidth + "px";
            let sliderStart = dateStart.offsetLeft + dateStart.offsetWidth;
            slider.style.left = sliderStart + "px";
            dateStop.style.left = (sliderStart + sliderWidth + canvasDiv.offsetLeft) + "px";
            slider.style.width = sliderWidth + "px";
            let start = sliderWidth * m_sgraph.m_Start / m_sgraph.m_nData + sliderStart;
            let stop = sliderWidth * m_sgraph.m_Stop / m_sgraph.m_nData + sliderStart - barStop.offsetWidth;
            barStart.style.left = start + "px";
            barStop.style.left = stop + "px";
            bar.style.left = (start + barStart.offsetWidth - 10) + "px";
            bar.style.width = (stop + 10 - (start + barStart.offsetWidth - 10)) + "px";
            start = Math.floor(m_sgraph.m_Start);
            if (start >= m_sgraph.m_nData)
                start = m_sgraph.m_nData - 1;
            if (start >= 0)
                dateStart.innerText = m_sgraph.m_DateTime[start];
            stop = Math.floor(m_sgraph.m_Stop);
            if (stop >= m_sgraph.m_nData)
                stop = m_sgraph.m_nData - 1;
            if (stop >= 0)
                dateStop.innerText = m_sgraph.m_DateTime[stop];
            m_profit = parseFloat(m_sgraph.m_profit).toFixed(2);
            bar.innerText = (m_profit >= 0) ? "Profit:$" + m_profit : "Profit:-$" + (-m_profit);
        }

        let info = document.getElementById("sg_info");
        if (info) {
            if (document.getElementById("param_cb34").checked) {
                info.innerHTML = m_sgraph.getPerfomance();
            } else {
                info.innerHTML = "";
            }
        }
    }
    function setGraphStartStop(dStart, dStop) {
        if (m_sgraph.m_nData > 1 && m_sgraph && m_sgraph.m_Stop > m_sgraph.m_Start) {
            m_sgraph.m_Start = dStart * m_sgraph.m_nData;
            m_sgraph.m_Stop = dStop * m_sgraph.m_nData;
            m_sgraph.myDrawImage();
        }
    }
    function hideParamsWindow(pid) {
        if (isNaN(pid)) {
            if (!pid.startsWith("sg_input"))
                return;
            ind = pid.replace(/\D/g, "");
        }
        else
            ind = pid;
        let params = getParamsWindowByInd(ind);
        /*
        if (parseInt(params.style.zIndex) != m_zIndexParamsMax && elem.endsWith("_div"))
            paramsWindowToTop(ind);
        else {
            let op = params.style.opacity;
            params.style.opacity = (!op) ? 0.6 : (op > 0.4) ? op - 0.4 : 1.0;
        }
        */
        params.style.display = "none";
        m_showParamWindow[ind - 1] = false;
        saveParamsWindowPos(ind)
    }
    function mouseOnGraphClick(e) {
        let clickX = getMouseX(e);
        let clickY = getMouseY(e);
        if (m_sgraph) {
            let n = m_sgraph.m_nData;
            let canv = document.getElementById('sg_canvas');
            if (canv && m_sgraph.m_nData > 1 && m_sgraph.m_optimize == 0) {
                let nCorr = m_CorrLength < 1 ? 2 : m_CorrLength * 2;
                let d = m_sgraph.m_Stop - m_sgraph.m_Start - 1;
                let iw = canv.offsetWidth;
                let x = clickX - canv.offsetLeft;
                let ind0 = Math.floor(m_sgraph.m_Start + (x * d) / (iw - 1) + 0.5);
                let y = clickY - canv.offsetTop;
                let ih = canv.offsetHeight;
                let ind1 = Math.floor(ind0 - nCorr - (ih - y - 1) * (MATRIX_Y - 1) / (ih - 1));
                if (m_sgraph.m_Index0 != ind0 || m_sgraph.m_Index1 != ind1) {
                    m_sgraph.m_Index0 = ind0;
                    m_sgraph.m_Index1 = ind1;

                    // Select corresponding row in historical data table
                    let tableFrame = document.getElementById("historical_data");
                    let table = document.getElementById("historical_data_table");
                    if (tableFrame && table && table.rows && table.rows.length > 0 && m_sgraph.m_Index0 >= 0 && m_sgraph.m_Index0 < table.rows.length - 1) {
                        let trows = table.rows;
                        for (let i = 0; i < trows.length; i++) {
                            let row = table.rows[i];
                            if (row.style != "")
                                row.style = "";
                        }
                        let row = table.rows[m_sgraph.m_Index0 + 1];
                        row.style = "background-color: #a0ffa0";
                        tableFrame.scrollTop = row.offsetTop - (tableFrame.clientHeight / 2);
                        row = row;
                    }
                    m_sgraph.myDrawImage();
                }
            }
        }
    }
    function updateSlidersValues() {
        let vrange = document.getElementById("param_sl2"); // m_OrderScale
        let vlabel = document.getElementById("param_sla2"); // m_OrderScale
        vrange.value = m_OrderScale.toString();
        vlabel.innerText = vrange.value;
        vrange = document.getElementById("param_sl4"); // m_MaxLost
        vlabel = document.getElementById("param_sla4"); // m_MaxLost
        vrange.value = m_MaxLost.toString();
        vlabel.innerText = vrange.value;
        vrange = document.getElementById("param_sl3"); // m_AverageLength
        vlabel = document.getElementById("param_sla3"); // m_AverageLength
        vrange.value = m_AverageLength.toString();
        vlabel.innerText = vrange.value;
        vrange = document.getElementById("param_sl1"); // m_CorrLength
        vlabel = document.getElementById("param_sla1"); // m_CorrLength
        vrange.value = m_CorrLength.toString();
        vlabel.innerText = vrange.value;
    }
    function updateSlidersLabels() {
        let vlabel = document.getElementById("param_sla2"); // m_OrderScale
        vlabel.innerText = m_OrderScale.toString();
        vlabel = document.getElementById("param_sla4"); // m_MaxLost
        vlabel.innerText = m_MaxLost.toString();
        vlabel = document.getElementById("param_sla3"); // m_AverageLength
        vlabel.innerText = m_AverageLength.toString();
        vlabel = document.getElementById("param_sla1"); // m_CorrLength
        vlabel.innerText = m_CorrLength.toString();
        vlabel = document.getElementById("param_sla5"); // m_nDaysToKeepOpen
        vlabel.innerText = m_nDaysToKeepOpen.toString();
        vlabel = document.getElementById("param_sla6"); // m_Spread
        vlabel.innerText = m_Spread.toString();
    }
    function onOptimizationEnd() {
        updateSlidersValues();
        saveParamsValues();
        document.getElementById("sg_optimize").innerText = "Optimize";
        m_sgraph.processLoadedData();
    }
    function getOptimizeRanges() {
        let ind = 0;
        let vrange = document.getElementById("param_sl2"); // m_OrderScale
        m_sgraph.m_ranges[ind + 0] = parseFloat(vrange.min);
        m_sgraph.m_ranges[ind + 1] = parseFloat(vrange.max);
        m_sgraph.m_ranges[ind + 2] = parseFloat(vrange.step);
        ind += 3;
        vrange = document.getElementById("param_sl4"); // m_MaxLost
        m_sgraph.m_ranges[ind + 0] = parseFloat(vrange.min);
        m_sgraph.m_ranges[ind + 1] = parseFloat(vrange.max);
        m_sgraph.m_ranges[ind + 2] = parseFloat(vrange.step);
        ind += 3;
        vrange = document.getElementById("param_sl3"); // m_AverageLength
        m_sgraph.m_ranges[ind + 0] = parseFloat(vrange.min);
        m_sgraph.m_ranges[ind + 1] = parseFloat(vrange.max);
        m_sgraph.m_ranges[ind + 2] = parseFloat(vrange.step);
        if (m_bKnownOrder || m_bFixedOrder) {
            m_sgraph.m_ranges[ind + 0] = m_AverageLength;
            m_sgraph.m_ranges[ind + 1] = m_sgraph.m_ranges[ind + 0];
        }
        ind += 3;
        vrange = document.getElementById("param_sl1"); // m_CorrLength
        m_sgraph.m_ranges[ind + 0] = parseFloat(vrange.min);
        m_sgraph.m_ranges[ind + 1] = parseFloat(vrange.max);
        m_sgraph.m_ranges[ind + 2] = parseFloat(vrange.step);
        if (m_bKnownOrder || m_bFixedOrder) {
            m_sgraph.m_ranges[ind + 0] = m_CorrLength;
            m_sgraph.m_ranges[ind + 1] = m_sgraph.m_ranges[ind + 0];
        }
    }
    function onClickedOptimize() {
        if (m_sgraph && m_sgraph.m_nData > 1) {
            if (m_sgraph.startOrStopOptimization())
                document.getElementById("sg_optimize").innerText = "Stop Optimization";
            else
                document.getElementById("sg_optimize").innerText = "Optimize";
        }
    }
    function paramSliderChange(indDiv, indRange) {
        let tableId = getParamsWindowByInd(indDiv);
        let labelId = document.getElementById("param_sla" + indRange);
        let rangeId = document.getElementById("param_sl" + indRange);
        if (parseInt(tableId.style.zIndex) != m_zIndexParamsMax) {
            paramsWindowToTop(tableId.id);
        }
        labelId.innerText = rangeId.value;
        showNewSGraph(); // fast enough to run for each change, no "setTimeout" needed
    }
    function getDateIndex(dateStr) {
        let ind = -1;
        let dateStr3 = dateStr.split(/\/|-/);
        let dateStr0 = parseInt(dateStr3[0]);
        let dateStr1 = parseInt(dateStr3[1]);
        let dateStr2 = parseInt(dateStr3[2]);
        if (!dateStr3 || dateStr3.length != 3)
            return ind;
        for (let i = 0; i < m_sgraph.m_nData; i++) {
            let date3 = m_sgraph.m_DateTime[i].split(/\/|-/);
            if (date3.length == 3) {
                let date0 = parseInt(date3[0]);
                let date1 = parseInt(date3[1]);
                let date2 = parseInt(date3[2]);
                if ((date0 === dateStr0 && date1 === dateStr1 && date2 === dateStr2) ||
                    (date0 === dateStr2 && date1 === dateStr0 && date2 === dateStr1))
                    return i;
            }
        }
        return ind;
    }
    function resutlTableRowClicked(elm) {
        let tm = new Date();
        let timeSinceastClick = tm.getTime() - m_lastArchiveClickTime.getTime();
        m_lastArchiveClickTime = tm;
        if (timeSinceastClick > 300) // wait for double click
            return;
        let index = elm.rowIndex;
        let table = document.getElementById("sg_params_table");
        let fileName = document.getElementById("sg_filename");
        if (table && fileName) {
            let rowStr = getTableRowStr(table, index);
            let ind1 = 0;
            let ind2 = 0;
            let values = rowStr.split(",");
            if (m_numberOfColumnsInResult != values.length) {
                alert("Number of values " + values.length + " different than expected(" + m_numberOfColumnsInResult + ")");
            }
            else {
                if (values[1] !== fileName.innerText || m_sgraph.m_nData < 1) {
                    alert("Please load file " + values[1] + " first");
                }
                else {
                    ind1 = getDateIndex(values[2]);
                    if (ind1 < 0) {
                        alert("Can't find date " + values[2] + " in file " + values[1] + ", SW will use the first date");
                        ind1 = 0;
                    }
                    ind2 = getDateIndex(values[3]);
                    if (ind2 < 0) {
                        alert("Can't find date " + values[3] + " in file " + values[1] + ", SW will use the last date");
                        ind2 = m_sgraph.m_nData - 1;
                    }
                    if (ind2 < ind1) {
                        alert("Incorrect Start/Stop dates, SW will use the first and last dates from file " + values[1]);
                        ind1 = 0;
                        ind2 = m_sgraph.m_nData - 1;
                    }
                    m_sgraph.m_Start = ind1;
                    m_sgraph.m_Stop = ind2;
                    for (let i = 1; i <= m_numberOfColumnsInResult - 6; i++) {
                        let vrange = document.getElementById("param_sl" + i);
                        let vlabel = document.getElementById("param_sla" + i);
                        vrange.value = values[i + 3];
                        vlabel.innerText = vrange.value;
                    }
                    let mode = values[m_numberOfColumnsInResult - 2];
                    m_bFixedOrder = (mode == "Fixed order") ? true : false;
                    m_bKnownOrder = (mode == "Known price") ? true : false;

                    let trows = table.rows;
                    for (let i = 0; i < trows.length; i++) {
                        let row = table.rows[i];
                        if (row.style != "")
                            row.style = "";
                    }
                    let row = table.rows[index];
                    row.style = "background-color: #a0ffa0";

                    showNewSGraph();
                }
            }
        }
    }
    function storeTableIfNeeded() {
        if (m_archiveTableChangedIndex > 0) {
            saveParamsTable(m_archiveTableChangedIndex);
            m_archiveTableChangedIndex = -1;
        }
    }
    function tableChanged(elm) {
        if (elm && elm.tabIndex) {
            m_archiveTableChangedIndex = elm.tabIndex;
        }
    }
    /* When the user clicks on the button, 
    toggle between hiding and showing the dropdown content */
    function showHideShowMenu() {
        document.getElementById("myDropdown").classList.toggle("show");
    }
    function getParamsWindowIdByInd(ind) {
        return "sg_input_params" + ind + "_div";
    }
    function getParamsWindowByInd(ind) {
        return document.getElementById(getParamsWindowIdByInd(ind));
    }
    function processWindowOnClick(event) {
        // Close the dropdown if click is outside
        if (event.target.id != "hide_show_params") {
            let dropdowns = document.getElementsByClassName("dropdown-content");
            for (let i = 0; i < dropdowns.length; i++) {
                let openDropdown = dropdowns[i];
                if (openDropdown.classList.contains('show')) {
                    openDropdown.classList.remove('show');
                }
            }
        }
        if (event.target.name === "sg_pd_menu") {
            let ind = m_params_windows.indexOf(event.target.innerText);
            if (ind == 0)
                showHideAllParamsWindows(true);
            else if (ind == m_nParamsWindows + 1)
                showHideAllParamsWindows(false);
            else if (ind > 0 && m_sgraph && m_sgraph.m_nData > 0) {
                if (m_showParamWindow[ind - 1])
                    hideParamsWindow(ind);
                else {
                    m_showParamWindow[ind - 1] = true;
                    paramsWindowToTop(ind);
                }
            }
        }
    }
    return {
        onZoomIn: onZoomIn,
        onZoomOut: onZoomOut,
        handleFile: handleFile,
        showSGraph: showSGraph,
        onLoadSetup: onLoadSetup,
        resizeGraph: resizeGraph,
        tableChanged: tableChanged,
        showNewSGraph: showNewSGraph,
        copyArhiveTable: copyArhiveTable,
        clickArhiveTable: clickArhiveTable,
        pasteArhiveTable: pasteArhiveTable,
        showHideShowMenu: showHideShowMenu,
        hideParamsWindow: hideParamsWindow,
        copyResultToTable: copyResultToTable,
        rejectArhiveInput: rejectArhiveInput,
        paramsWindowClick: paramsWindowClick,
        onClickedOptimize: onClickedOptimize,
        getOptimizeRanges: getOptimizeRanges,
        paramsWindowToTop: paramsWindowToTop,
        paramSliderChange: paramSliderChange,
        storeTableIfNeeded: storeTableIfNeeded,
        updateSlidersValues: updateSlidersValues,
        copyTableToClipboard: copyTableToClipboard,
        resutlTableRowClicked: resutlTableRowClicked,
        get m_cnv() { return m_cnv; }, set m_cnv(v) { m_cnv = v; },
        get m_ctx() { return m_ctx; }, set m_ctx(v) { m_ctx = v; },
        get m_Spread() { return m_Spread; }, set m_Spread(v) { m_Spread = v; },
        get MATRIX_Y() { return MATRIX_Y; }, set MATRIX_Y(v) { MATRIX_Y = v; },
        get m_bShowOpen() { return m_bShowOpen; }, set m_bShowOpen(v) { m_bShowOpen = v; },
        get m_bShowCorr() { return m_bShowCorr; }, set m_bShowCorr(v) { m_bShowCorr = v; },
        get m_bShowOrder() { return m_bShowOrder; }, set m_bShowOrder(v) { m_bShowOrder = v; },
        get m_bKnownOrder() { return m_bKnownOrder; }, set m_bKnownOrder(v) { m_bKnownOrder = v; },
        get m_bShowVolume() { return m_bShowVolume; }, set m_bShowVolume(v) { m_bShowVolume = v; },
        get m_bShowProfit() { return m_bShowProfit; }, set m_bShowProfit(v) { m_bShowProfit = v; },
        get m_bFixedOrder() { return m_bFixedOrder; }, set m_bFixedOrder(v) { m_bFixedOrder = v; },
        get m_bShowCorrPos() { return m_bShowCorrPos; }, set m_bShowCorrPos(v) { m_bShowCorrPos = v; },
        get m_bShowAverage() { return m_bShowAverage; }, set m_bShowAverage(v) { m_bShowAverage = v; },
        get m_MaxLost() { return m_MaxLost; }, set m_MaxLost(v) { m_MaxLost = v; },
        get m_nDaysToKeepOpen() { return m_nDaysToKeepOpen; }, set m_nDaysToKeepOpen(v) { m_nDaysToKeepOpen = v; },
        get m_bShowPerfomance() { return m_bShowPerfomance; }, set m_bShowPerfomance(v) { m_bShowPerfomance = v; },
        get m_bCalculatedOrder() { return m_bCalculatedOrder; }, set m_bCalculatedOrder(v) { m_bCalculatedOrder = v; },
        get m_bShowCandleStick() { return m_bShowCandleStick; }, set m_bShowCandleStick(v) { m_bShowCandleStick = v; },
        get m_CorrLength() { return m_CorrLength; }, set m_CorrLength(v) { m_CorrLength = v; },
        get m_OrderScale() { return m_OrderScale; }, set m_OrderScale(v) { m_OrderScale = v; },
        get m_AverageLength() { return m_AverageLength; }, set m_AverageLength(v) { m_AverageLength = v; },
        get m_AverageLengthMax() { return m_AverageLengthMax; }, set m_AverageLengthMax(v) { m_AverageLengthMax = v; },
        get m_bHighlightCorrIntervals() { return m_bHighlightCorrIntervals; }, set m_bHighlightCorrIntervals(v) { m_bHighlightCorrIntervals = v; }
    };
}
