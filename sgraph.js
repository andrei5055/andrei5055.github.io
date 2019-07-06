// SGraph - Create bitmap with correlation colormap and optional graphics
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
function sgraph()
{
    /* following "enums" must be on top of all other attributes */
    var ind = 0;
    var INDEX_OPEN = ind; ind++;
    var INDEX_HIGH = ind; ind++;
    var INDEX_LOW = ind; ind++;
    var INDEX_CLOSE = ind; ind++;
    var INDEX_ADJ_CLOSE = ind; ind++;
    var INDEX_VOLUME = ind; ind++;
    var INDEX_AVERAGE = ind; ind++;
    var INDEX_PROFIT = ind; ind++;
    var INDEX_ORDER = ind; ind++;
    var INDEX_MAX_CORR = ind; ind++;
    var INDEX_MAX_CORR_DT = ind; ind++;
    var NDATA = ind;

    ind = 0;
    var INDEX_OPEN_FP = ind; ind++;
    var INDEX_OPEN_SPREAD_FP = ind; ind++;
    var INDEX_HIGH_FP = ind; ind++;
    var INDEX_HIGH_SPREAD_FP = ind; ind++;
    var INDEX_LOW_FP = ind; ind++;
    var INDEX_LOW_SPREAD_FP = ind; ind++;
    var INDEX_CLOSE_FP = ind; ind++;
    var INDEX_CLOSE_SPREAD_FP = ind; ind++;
    var INDEX_MAX_LOST_FP = ind; ind++;
    var INDEX_DELTA_FP = ind; ind++;
    var INDEX_DELTA2_FP = ind; ind++;
    var NDATA_FP = ind;

    ind = 0;
    var eNoOptimization = ind; ind++;
    var eStart1stOptimizationLoop = ind; ind++;
    var eRun1stOptimizationLoop = ind; ind++;
    var eStop1stOptimizationLoop = ind; ind++;
    var eStart2ndOptimizationLoop = ind; ind++;
    var eRun2ndOptimizationLoop = ind; ind++;
    var eStop2ndOptimizationLoop = ind; ind++;
    var eEndOfOptimization = ind; ind++;
    /* add new attributes after "enums" above */

    var HIGHLIGHT = 35;
    var MIN_LENGTH = 100;
    var N_C = 10;
    var dt1 = new Date();
    var dt2 = new Date();
    var dt3 = new Date();
    var dt4 = new Date();
    var dt5 = new Date();
    var dt6 = new Date();
    var dt7 = new Date();
    var dt8 = new Date();

    var onOptimizeEnd = null;
    var onDrawEnd = null;
    var m_image = null;
    var m_imageCanvas = null;
    var m_imWidth = 0;
    var m_imHeight = 0;
    var m_corrWidth = 0;
    var m_corrHeight = 0;
    var m_Index0 = -1;
    var m_Index1 = -1;

    var m_Start = 0;
    var m_Stop = 0;
    var m_nData = null;
    var m_ranges = new Array(0);
    var m_optimumProfit = 0.0;
    var m_optimize = eNoOptimization;

    var m_DateTime = new Array();

    var progress = 0;

    var m_xMouseDown = 0;
    var m_yMouseDown = 0;
    var m_data = new Array();
    var m_dataFP = new Array();
    var m_dataCorr = new Array();
    var m_Corr = new Array();
    var m_gName = new Array();
    var m_gMin = new Array();
    var m_gMax = new Array();
    var m_gDraw = new Array();
    var m_gColors = new Array();

    var m_dMin = 0.0;
    var m_dMax = 0.0;
    var m_dMinPr = 0.0;
    var m_dMaxPr = 0.0;

    var m_prevMaxCorr = 0.0;
    var m_prevMaxCorrPos = -1;
    var m_optimizeTimeout = 0;
    var m_optimum_PredScale = 0;
    var m_optimum_MaxLost = 0;
    var m_optimum_ave = 0;
    var m_optimum_ncorr = 0;
    var m_profit = 0.0;
    zoom = function(iZoom)
    {
        if (m_nData > 0) {
            let minLng = m_nData > MIN_LENGTH ? MIN_LENGTH : m_nData;
            if (m_Stop > m_nData)
                m_Stop = m_nData;
            if (m_Start < 0)
                m_Start = 0;
            let n = Math.floor(m_Stop - m_Start);
            if (n > m_nData)
                n = m_nData;

            let n1 = iZoom >= 0 ? n * 2 / 3 : n * 3 / 2;

            if (n1 < minLng)
                n1 = minLng;
            if (n1 > m_nData)
                n1 = m_nData;

            if (m_Index0 >= m_Start && m_Index0 <= m_Stop) {
                m_Start = m_Index0 - (m_Index0 - m_Start) * n1 / n;
                m_Stop = m_Index0 + (m_Stop - m_Index0) * n1 / n;
            }
            else
                m_Start = m_Stop - n1;
            if (m_Start < 0) {
                m_Start = 0;
                m_Stop = m_Start + n1;
            }
        }
    }
    drawCandlestick = function(ix1, iww, open, high, low, close)
    {
        let closeOpenMax = Math.max(close, open);
        let closeOpenMin = Math.min(close, open);
        let iy = Math.floor(high);
        let iy1 = Math.floor(closeOpenMax);
        let iww2 = Math.floor(iww / 2);
            
        let urRect = 0;
        let ugRect = 255;
        let ubRect = 0;
        if (close < open)
        {
            urRect = 240;
            ugRect = 100;
            ubRect = 50;
        }
        sg_main.m_ctx.beginPath();
        sg_main.m_ctx.moveTo(ix1, m_imHeight - 1 - iy);
        sg_main.m_ctx.lineTo(ix1, m_imHeight - 1 - iy1);
        sg_main.m_ctx.stroke();

        iy = iy1;
        iy1 = Math.floor(closeOpenMin);
        sg_main.m_ctx.globalAlpha = 0.7;
        sg_main.m_ctx.fillStyle = RGB2Hex(urRect, ugRect, ubRect);
        sg_main.m_ctx.fillRect(ix1 - iww2, m_imHeight - 1 - iy, iww, iy - iy1 + 1);
        sg_main.m_ctx.globalAlpha = 1.0;

        sg_main.m_ctx.beginPath();
        iy = iy1;
        iy1 = Math.floor(low);
        sg_main.m_ctx.moveTo(ix1, m_imHeight - 1 - iy);
        sg_main.m_ctx.lineTo(ix1, m_imHeight - 1 - iy1);
        sg_main.m_ctx.stroke();
    }
    calculateSumXY = function(iShift)
    {
        // Calculate Sum(X*Y)
        let s2 = 0.0;
        let x0 = 0.0, x1;
        let y0 = 0.0, y1;
        let indx = (iShift+1) * NDATA;
        let indy = NDATA;
        let j = (iShift+1)*N_C+2;
        for (let i = iShift+1; i < m_nData; i++, j += N_C, indx += NDATA, indy += NDATA)
        {
            x1 = m_data[indx + INDEX_HIGH];
            y1 = m_data[indy + INDEX_HIGH];
            x1 -= m_data[indx + INDEX_AVERAGE];
            y1 -= m_data[indy + INDEX_AVERAGE];
            s2 += x1*y1;
            m_dataCorr[j] = s2;
        }
    }
    calculateCorrelation = function(j, nCorr)
    {
        let iCorr = j * m_nData;
        let iShift = j+nCorr;
    
        calculateSumXY(iShift);

        let indCorr = iShift+nCorr+1;
        if (indCorr < m_nData)
        {
            for (let i = 0; i < indCorr; i++)
                m_Corr[iCorr+i] = 0.0;
        }
        else
            iShift = iShift;

        let indx = indCorr * N_C;
        let indy = indx - iShift * N_C;
        let indxs = indx - nCorr * N_C;
        let indys = indy - nCorr * N_C;
        for (let i = indCorr; i < m_nData; i++, indx += N_C, indy += N_C, indxs += N_C, indys += N_C)
        {
            let a = nCorr *(m_dataCorr[indx+2] - m_dataCorr[indxs+2]);
            let bx = m_dataCorr[indx] - m_dataCorr[indxs];
            let by = m_dataCorr[indy] - m_dataCorr[indys];
            let b = bx * by;
            let cxx = nCorr *(m_dataCorr[indx+1] - m_dataCorr[indxs+1]);
            let cx = bx * bx ;
            let cc = cxx - cx;
            let dyy = nCorr *(m_dataCorr[indy+1] - m_dataCorr[indys+1]);
            let dy = by * by ;
            let dd = dyy - dy;
            let ee = cc * dd;
            let corr = (ee > 0.0000000001) ? (a - b) / Math.sqrt(ee) : 0.0;
            m_Corr[iCorr+i] = corr;
            if (m_minCorr > corr)
                m_minCorr = corr;
            if (m_maxCorr < corr)
                m_maxCorr = corr;
        }
    }
    getInterpolatedValueInt = function (ind, width, start, stop) {
        return (Math.floor((ind - start) * (width - 1) / (stop - start - 1)));
    }
    getInterpolatedValue = function (ind, width, start, stop) {
        return ((ind - start) * (width - 1) / (stop - start - 1));
    }
    orderString = function()
    {
        let orderStr = "Not available";
        let ind = m_nData;
        if (ind > 0) {
            if (sg_main.m_bKnownOrder)
                orderStr = "'Known tomorrow price' order for next day (after " + m_DateTime[ind - 1] + "):\nSale limit = tomorrow High (-:)";
            else {
                let indy = m_data[(ind - 1) * NDATA + INDEX_MAX_CORR_DT];
                if (indy >= 0 && indy < sg_main.MATRIX_Y) {
                    let corrInd = m_nData * indy + ind - 1;
                    let nCorr = sg_main.m_CorrLength * 2;
                    let corr = m_Corr[corrInd];
                    let dOrder = getOrder(corr, ind, nCorr, indy);
                    let openScale = 1.0 + sg_main.m_OrderScale * dOrder;

                    orderStr = sg_main.m_bFixedOrder ? "Fixed order " : "Calculated order ";
                    orderStr += "for next day (after " + m_DateTime[m_nData - 1] + "):\n";
                    orderStr += "Sale limit = Open x " + parseFloat(openScale).toFixed(3);
                }
            }
        }
        return orderStr;
    }
    getProfit = function()
    {
        let s = 0.0;
        let dMin = 999999.0;
        let dMax =-999999.0;
        let dMinPr = 1e20;
        let dMaxPr =-1e20;
        let spread = sg_main.m_Spread;
        let buyPrice = 0.0;
        let bOrder = false;
        let sellLimit = buyPrice;
        let prof = 0.0;
        let nOrders = 0;
        let close = 0.0;
        let todayProfit = 0.0;
        let nDaysOpenPos = 0;
        let nDaysToKeepOpen = sg_main.m_nDaysToKeepOpen;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        if (iStart < 2)
            iStart = 2;
        if (iStop < m_nData)
            iStop++;
        for (let i = iStart, indData = iStart * NDATA_FP; i < iStop; i++, indData += NDATA_FP)
        {
            let open = m_dataFP[indData + INDEX_OPEN_FP];
            let highSpread = m_dataFP[indData + INDEX_HIGH_SPREAD_FP];
            close = m_dataFP[indData + INDEX_CLOSE_FP];
            let dMaxLost = m_dataFP[indData + INDEX_MAX_LOST_FP];
            let delta = m_dataFP[indData + INDEX_DELTA_FP];
  
            todayProfit = 0.0;
        
            if (!bOrder)
            {
                if (delta > spread)
                {
                    buyPrice = open + spread;
                    sellLimit = open + delta;
                    nOrders++;
                    bOrder = true;
                    nDaysOpenPos = 0;
                }
            }

            if (bOrder)
            {
                if (dMaxLost != 0.0)
                {
                    let openSpread = m_dataFP[indData + INDEX_OPEN_SPREAD_FP];
                    if (buyPrice - openSpread > dMaxLost)
                    {
                        todayProfit = openSpread - buyPrice;
                        bOrder = false;
                    }
                    else
                    {
                        if (sellLimit <= highSpread && open > close)
                        {
                            // Assume (if open > close) sell limit triggered before lost limit 
                        }
                        else if(buyPrice - m_dataFP[indData + INDEX_LOW_SPREAD_FP] > dMaxLost)
                        {
                            todayProfit = -dMaxLost;
                            bOrder = false;
                        }
                    }
                }
                if (bOrder)
                {
                    if (sellLimit <= highSpread)
                    {
                        todayProfit = sellLimit - buyPrice;
                        bOrder = false;
                    }
                    else
                    {
                        if (nDaysOpenPos < nDaysToKeepOpen)
                        {
                            nDaysOpenPos++;
                        }
                        else
                        {
                            prof = m_dataFP[indData + INDEX_CLOSE_SPREAD_FP] - buyPrice;
                            // close position at the end of the day
                            todayProfit = prof;
                            bOrder = false;
                        }
                    }
                }
            }
            s += todayProfit;
            m_data[indData + INDEX_PROFIT] = s;
            if (dMin > delta)
                dMin = delta;
            if (dMax < delta)
                dMax = delta;
            if (dMinPr > s)
                dMinPr = s;
            if (dMaxPr < s)
                dMaxPr = s;
        }
        if (iStop == m_nData && bOrder)
        {
            todayProfit = close - spread - buyPrice;
            s += todayProfit;
        }
        m_dMin = dMin;
        m_dMax = dMax;
        m_dMinPr = dMinPr;
        m_dMaxPr = dMaxPr;
        m_profit = s;

        return s;
    }
    RGB2Hex = function(r, g, b) {
        r = r.toString(16);
        g = g.toString(16);
        b = b.toString(16);

        if (r.length == 1)
            r = "0" + r;
        if (g.length == 1)
            g = "0" + g;
        if (b.length == 1)
            b = "0" + b;

        return "#" + r + g + b;
    }
    drawGrpah = function(name, ind, r, g, b, dMin, dMax, bDraw)
    {
        m_gName[ind] = name;
        m_gMin[ind] = dMin;
        m_gMax[ind] = dMax;
        m_gDraw[ind] = bDraw;
        m_gColors[ind] = RGB2Hex(r, g, b);
        if (!bDraw)
            return;
        let ix1;
        let pt1x, pt2x, pt1y, pt2y;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        if (iStop < m_Stop)
            iStop++;
        if (iStop > m_nData)
            iStop = m_nData;
        if (iStop - iStart < 2)
            return;
    
        let iww = Math.floor((m_imWidth - 1) / (iStop - iStart - 1));
        if (iww > 5)
            iww = Math.floor(iww/3)*2+1;
        iww = Math.floor((iww - 1) / 2) * 2 + 1;

        let bShowCandlestick = false;
        if (Math.abs(dMax - dMin) < 1e-10)
        {
            dMax += 1.0;
            dMin -= 1.0;
        }
        let dScale = dMax - dMin;
        if (dScale <= 0.0)
            dScale = 1.0;
        dScale = (m_imHeight-1) / dScale;
        pt2x = getInterpolatedValue(iStart, m_imWidth, m_Start, m_Stop);
        pt2y = (m_data[iStart*NDATA + ind] - dMin) * dScale;

        if (ind == INDEX_OPEN)
            bShowCandlestick = sg_main.m_bShowCandleStick;
        let dataScale = 1.0;
        if (ind == INDEX_ORDER)
            dataScale = sg_main.m_OrderScale;

        sg_main.m_ctx.strokeStyle = m_gColors[ind];

        if (bShowCandlestick) {
            sg_main.m_ctx.lineWidth = 1;
            for (let i = iStart; i < iStop; i++)
            {
                let j = i*NDATA;
                ix1 = getInterpolatedValueInt(i, m_imWidth, m_Start, m_Stop);
                let open  = m_data[j + INDEX_OPEN];
                let high  = m_data[j + INDEX_HIGH];
                let low   = m_data[j + INDEX_LOW];
                let close = m_data[j + INDEX_CLOSE];
                open = ((open - dMin) * dScale);
                high = ((high - dMin) * dScale);
                low = ((low - dMin) * dScale);
                close = ((close - dMin) * dScale);
                drawCandlestick(ix1, iww, open, high, low, close);
            }
            if (!sg_main.m_bShowOpen)
                return;
        }
        sg_main.m_ctx.beginPath();
        sg_main.m_ctx.lineWidth = 3;
        for (let i = iStart; i < iStop; i++)
        {
            let j = i*NDATA;
            let v;
            if (ind == INDEX_ORDER)
                v = m_data[j + ind] * dataScale * m_data[j + INDEX_HIGH];
            else
                v = (m_data[j + ind]) * dataScale;
            
            pt1x = pt2x;
            pt1y = pt2y;
            pt2x = getInterpolatedValue(i, m_imWidth, m_Start, m_Stop);
            pt2y = Math.floor((v - dMin) * dScale);
            if (i == iStart)
                sg_main.m_ctx.moveTo(pt2x, m_imHeight - 1 - pt2y);
            else
                sg_main.m_ctx.lineTo(pt2x, m_imHeight - 1 - pt2y);
        }
        sg_main.m_ctx.stroke();
    }
  
    getVolumeMinMax = function()
    {
        let dMin = 1e20;
        let dMax = -1e20;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        for (let i = iStart, ii = iStart*NDATA; i < iStop; i++, ii += NDATA)
        {
            if (dMin > m_data[ii+INDEX_VOLUME])
                dMin = m_data[ii+INDEX_VOLUME];
            if (dMax < m_data[ii+INDEX_VOLUME])
                dMax = m_data[ii+INDEX_VOLUME];
        }
        return { dMin:dMin, dMax:dMax };
    }
    getStockMinMax = function()
    {
        let dMin = 1e20;
        let dMax = -1e20;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        for (let i = iStart, ii = iStart*NDATA; i < iStop; i++, ii += NDATA)
        {
            if (dMin > m_data[ii+INDEX_LOW])
                dMin = m_data[ii+INDEX_LOW];
            if (dMax < m_data[ii+INDEX_HIGH])
                dMax = m_data[ii+INDEX_HIGH];
        }
        return { dMin: dMin, dMax: dMax };
    }
    highlightCorrIntervals = function(nCorr)
    {
        let iw = m_imWidth;
        let ih = m_imHeight;
        if (sg_main.m_bHighlightCorrIntervals)
        {
            sg_main.m_ctx.fillStyle = "rgba(255, 255, 255, 0.3)";
            let ihx = [m_Index0 - nCorr, m_Index0, m_Index1 - nCorr, m_Index1];
            for (let i = 0; i < 4; i += 2)
            {
                if (ihx[i] >= 0 && ihx[i+1] >= ihx[i])
                {
                    let ix0 = getInterpolatedValueInt(ihx[i], iw, m_Start, m_Stop);
                    let ix1 = getInterpolatedValueInt(ihx[i + 1], iw, m_Start, m_Stop);
                    if (ix0 < 0)
                        ix0 = 0;
                    if (ix1 >= iw)
                        ix1 = iw - 1;
                    if (ix0 <= ix1) {
                        sg_main.m_ctx.fillRect(ix0, 0, ix1 - ix0 + 1, ih);
                    }
                }
            }
        }
        sg_main.m_ctx.lineWidth = 1;
        sg_main.m_ctx.globalAlpha = 0.3;
        sg_main.m_ctx.strokeStyle = "white";
        let ix1 = getInterpolatedValueInt(m_Index0, iw, m_Start, m_Stop);
        let iy1 = getInterpolatedValueInt(m_Index0 - m_Index1 - nCorr, ih, 0, sg_main.MATRIX_Y) - 1;

        sg_main.m_ctx.beginPath();
        if (ix1 >= 0 && ix1 < iw) {
            sg_main.m_ctx.moveTo(ix1, 0);
            sg_main.m_ctx.lineTo(ix1, ih);
            sg_main.m_ctx.stroke();
            sg_main.m_ctx.beginPath();
            sg_main.m_ctx.moveTo(ix1-5, ih - 1 - iy1);
            sg_main.m_ctx.lineTo(ix1 + 5, ih - 1 - iy1);
            sg_main.m_ctx.stroke();
        }
        sg_main.m_ctx.globalAlpha = 1.0;
    }
    scaleImage = function (image_data_src, width, height, w, h, dStart, dStop) {
        let imageCanvas_data = m_imageCanvas.data;
        let imgBytesPerPixel = 4;
        let iwInp = w * imgBytesPerPixel;
        let iwOut = width * imgBytesPerPixel;
        let ddx = (dStop - 1.0 - dStart) / (width - 1);
        let ddy = (h - 1) / (height - 1);
        let dy = 0.0;
        let dZoomOriginX = dStart;
        for (let j = 0; j < height; j++) {
            let dx = dZoomOriginX;
            let jj = Math.floor(dy);
            let yd = dy - jj;
            let cOut = (height - 1 - j) * iwOut;
            let cIn = jj * iwInp;

            if (jj >= h)
                break;
            if (jj == h - 1)
                iwInp = 0;
            for (let i = 0; i < width; i++) {
                let ii = Math.floor(dx);
                if (ii >= w)
                    break; // Temporary... Cannot happen
                let xd = dx - ii;
                let cI1 = cIn + ii * imgBytesPerPixel;
                let cI2 = cI1 + iwInp;
                if (ii == w - 1) {
                    imageCanvas_data[cOut] = (image_data_src[cI1] * (1.0 - yd) + image_data_src[cI2] * yd);
                    imageCanvas_data[cOut + 1] = (image_data_src[cI1 + 1] * (1.0 - yd) + image_data_src[cI2 + 1] * yd);
                    imageCanvas_data[cOut + 2] = (image_data_src[cI1 + 2] * (1.0 - yd) + image_data_src[cI2 + 2] * yd);
                    imageCanvas_data[cOut + 3] = 255;
                }
                else {
                    imageCanvas_data[cOut] = ((image_data_src[cI1] * (1.0 - xd) + image_data_src[cI1 + imgBytesPerPixel] * xd) * (1 - yd) + (image_data_src[cI2] * (1.0 - xd) + image_data_src[cI2 + imgBytesPerPixel] * xd) * yd);
                    imageCanvas_data[cOut + 1] = ((image_data_src[cI1 + 1] * (1.0 - xd) + image_data_src[cI1 + (imgBytesPerPixel + 1)] * xd) * (1 - yd) + (image_data_src[cI2 + 1] * (1.0 - xd) + image_data_src[cI2 + (imgBytesPerPixel + 1)] * xd) * yd);
                    imageCanvas_data[cOut + 2] = ((image_data_src[cI1 + 2] * (1.0 - xd) + image_data_src[cI1 + (imgBytesPerPixel + 2)] * xd) * (1 - yd) + (image_data_src[cI2 + 2] * (1.0 - xd) + image_data_src[cI2 + (imgBytesPerPixel + 2)] * xd) * yd);
                    imageCanvas_data[cOut + 3] = 255;
                }
                cOut += imgBytesPerPixel;
                dx += ddx;
            }
            dy += ddy;
        }
    }
    fillBitmap = function()
    {
        let iw = m_imWidth;
        let ih = m_imHeight;
        let ind = INDEX_OPEN;
        let r = 0;
        let g = 255;
        let b = 0;
        let ix, iy;
        let nCorr = sg_main.m_CorrLength * 2;
        if (nCorr < 2)
            nCorr = 2;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        if (iStop < m_Stop)
            iStop++;
        if (iStop > m_nData)
            iStop = m_nData;
        if (iStop - iStart < 2)
            return 0;
        let stockMinMax = getStockMinMax();
        let volumeMinMax = getVolumeMinMax();
        dMinA = stockMinMax.dMin;
        dMaxA = stockMinMax.dMax;
        let d = dMinA - dMinA;
        if (d <= 0.0)
            d = 1.0;
        d = (ih - 1) / d;
        setCorrelationMatrix();
        sg_main.m_ctx.clearRect(0, 0, m_imWidth, m_imHeight);
        scaleImage(m_image.data, iw, ih, m_image.width, m_image.height, m_Start, m_Stop);

        sg_main.m_ctx.putImageData(m_imageCanvas, 0, 0);

        drawGrpah("Open", ind, r, g, b, dMinA, dMaxA, sg_main.m_bShowOpen || sg_main.m_bShowCandleStick);
   
        highlightCorrIntervals(nCorr);
        
        getProfit();

        let dd = m_dMax - m_dMin;
        let dMin = m_dMin - dd * 0.5;
        let dMax = m_dMax + dd * 0.5;
        drawGrpah("Average",     INDEX_AVERAGE,     200, 200, 200, dMinA,  dMaxA,    sg_main.m_bShowAverage);
        drawGrpah("Volume",      INDEX_VOLUME,      200, 200, 200, volumeMinMax.dMin,  volumeMinMax.dMax,    sg_main.m_bShowVolume);
        drawGrpah("Order",       INDEX_ORDER,       127, 127, 255, dMin,   dMax,     sg_main.m_bShowOrder);
        drawGrpah("Profit",      INDEX_PROFIT,      255, 255, 255, m_dMinPr, m_dMaxPr,   sg_main.m_bShowProfit);
        drawGrpah("Correlation", INDEX_MAX_CORR,    255, 127, 127,   -1.0,   1.0,      sg_main.m_bShowCorr);
        drawGrpah("Maximum Correlation", INDEX_MAX_CORR_DT, 255, 0, 255, 0.0, sg_main.MATRIX_Y, sg_main.m_bShowCorrPos);

        return 1;
    }
    calculateOrderDelta = function ()
    {
        let dOrderScale = sg_main.m_OrderScale;
        for (let i = 0, indData = 0, outdData = 0; i < m_nData; i++, indData += NDATA, outdData += NDATA_FP)
        {
            let openRaw = m_data[indData + INDEX_OPEN];
            let order = m_data[indData + INDEX_ORDER] * openRaw;
            let delta = order * dOrderScale;
            m_dataFP[outdData + INDEX_DELTA_FP] = delta;
        }
    }
    calculateMaxLost = function()
    {
        let dMaxLostScale = sg_main.m_MaxLost * 0.01;
        for (let i = 0, indData = 0, outdData = 0; i < m_nData; i++, indData += NDATA, outdData += NDATA_FP)
        {
            m_dataFP[outdData + INDEX_MAX_LOST_FP]  = m_data[indData + INDEX_OPEN] * dMaxLostScale;
        }
    }
    getAdjustedBySpreadOHLC = function()
    {
        let spread = sg_main.m_Spread;
        for (let i = 0, indData = 0, outdData = 0; i < m_nData; i++, indData += NDATA, outdData += NDATA_FP)
        {
            let open = m_data[indData + INDEX_OPEN];
            let high = m_data[indData + INDEX_HIGH];
            let low = m_data[indData + INDEX_LOW];
            let close = m_data[indData + INDEX_CLOSE];
            let openSpread  = open  - spread;
            let highSpread  = high  - spread;
            let lowSpread   = low   - spread;
            let closeSpread = close - spread;
            m_dataFP[outdData + INDEX_OPEN_FP]         = open;
            m_dataFP[outdData + INDEX_HIGH_FP]         = high;
            m_dataFP[outdData + INDEX_LOW_FP]          = low;
            m_dataFP[outdData + INDEX_CLOSE_FP]        = close;
            m_dataFP[outdData + INDEX_OPEN_SPREAD_FP]  = openSpread;
            m_dataFP[outdData + INDEX_HIGH_SPREAD_FP]  = highSpread;
            m_dataFP[outdData + INDEX_LOW_SPREAD_FP]   = lowSpread;
            m_dataFP[outdData + INDEX_CLOSE_SPREAD_FP] = closeSpread;
        }
    }
    calculateArraysForCorrelation = function()
    {
        getAdjustedBySpreadOHLC();
        calculateMaxLost();
    }
    deleteCorrArrays = function()
    {
        m_Corr.length = null;
        m_dataCorr.length = null;
    }
    getCorrelationMatrix = function()
    {
        if(m_dataCorr.length == 0)
        {
            let length = m_nData * N_C;
            let nCorr = sg_main.m_CorrLength * 2;
            if (nCorr < 2)
                nCorr = 2;
            let s0 = 0.0;
            let s1 = 0.0;
            let indx = NDATA;
            let x;
            // Calculate Sum(X), and Sum(X*X)
            for (let i = 1, j = N_C; i < m_nData; i++, j += N_C, indx += NDATA)
            {
                x = m_data[indx + INDEX_HIGH] - m_data[indx + INDEX_AVERAGE];
                s0 += x;
                m_dataCorr[j] = s0;
                s1 += x*x;
                m_dataCorr[j+1] = s1;
            }

            m_minCorr = 1e+17;
            m_maxCorr = -1e+17;
            for (let j = 0; j < sg_main.MATRIX_Y; j++)
            {
                calculateCorrelation(j, nCorr);
            }
            if (Math.abs(m_minCorr-m_maxCorr) < 0.0001)
            {
                m_minCorr = -1.0;
                m_maxCorr = 1.0;
            }
        }
    }
    
    getAverage = function()
    {
        let dAve = 0;
        let nCorrAve = sg_main.m_AverageLength * 2 + 1;
    
        if (nCorrAve <= 0)
        {
            for (let i = 0, ii = 0; i < m_nData; i++, ii += NDATA)
                m_data[ii + INDEX_AVERAGE]  = 0.0;
        }
        else
        {
            dAve = m_data[INDEX_OPEN];
            let d3 = (nCorrAve+1.0) / (sg_main.m_AverageLengthMax * 2.0 + 2.0);
            d3 = 1.0 - d3;
            for (let i = 0, ind = 0; i < m_nData; i++, ind += NDATA)
            {
                let v = m_data[ind + INDEX_OPEN];
                let diff = Math.abs(v-dAve);
                if (dAve > 0.0)
                    diff = diff * d3 / dAve;
                dAve = diff < 1.0 ? dAve * (1.0 - diff) + v * diff : v;
                m_data[ind + INDEX_AVERAGE] = dAve;
            }
        }
    }
    allocateImage = function(nx, ny)
    {
        if (m_image == null || m_corrWidth != nx || m_corrHeight != ny) {
            m_corrWidth = nx;
            m_corrHeight = ny;
            return sg_main.m_ctx.getImageData(0, 0, nx, ny);
        }
        return m_image;
    }
    setCorrelationMatrix = function()
    {
        m_image = allocateImage(m_nData, sg_main.MATRIX_Y);
        if (m_Corr.length == 0)
        {
            let nCorr = sg_main.m_CorrLength * 2;
            if (nCorr < 2)
                nCorr = 2;

            getAverage();

            getCorrelationMatrix();

            m_prevMaxCorPos = -1;

            for (let i = 1; i < m_nData; i++)
            {
                calcCorrelationColormapColumn(m_nData - i, nCorr);
            }
            m_data[INDEX_ORDER] = m_data[NDATA + INDEX_ORDER];
            m_data[INDEX_ORDER + 1] = m_data[NDATA + INDEX_ORDER + 1];
            m_data[INDEX_ORDER + 2] = m_data[NDATA + INDEX_ORDER + 2];

            calculateOrderDelta();
        }
    }
    getOrder = function (corr, i0, nCorr, indyCorr)
    {
        let i1Start = i0 - nCorr;
        let iCorrMax = i1Start - indyCorr;
        let ind = sg_main.m_bKnownOrder ? i0 : iCorrMax;
        let indn = ind * NDATA;
        let dOrder = 0.0;
        if (iCorrMax >= 0 && ind >= 0 && ind < m_nData && i1Start - iCorrMax >= 0 && i1Start - iCorrMax < sg_main.MATRIX_Y) {
            if (sg_main.m_bFixedOrder) {
                dOrder = 0.05;
            }
            else {
                if (sg_main.m_bKnownOrder) {
                    dOrder = m_data[indn + INDEX_HIGH] - m_data[indn + INDEX_OPEN];
                    dOrder -= sg_main.m_Spread;
                    dOrder = Math.abs(m_data[indn + INDEX_OPEN]) > 0.0 ? dOrder / m_data[indn + INDEX_OPEN] : 0.0;
                }
                else {
                    indn += NDATA;
                    dOrder = m_data[indn + INDEX_HIGH] - m_data[indn + INDEX_OPEN];
                    dOrder -= sg_main.m_Spread;
                    if (dOrder > 0.0) {
                        dOrder = (m_data[indn + INDEX_OPEN] > 0.0) ? 1.0 / m_data[indn + INDEX_OPEN] : 0.0;
                        if (m_maxCorr > m_minCorr)
                            dOrder = dOrder * (corr - m_minCorr) / (m_maxCorr - m_minCorr);
                    }
                    else
                        dOrder = 0.0;
                }
            }
        }
        return dOrder;
    }
    calcCorrelationColormapColumn = function (i0, nCorr, dInitialCorrPos)
    {
        let nbpl = 4;
        let ih = m_image.height;
        let lw = m_image.width*nbpl;
        let index = i0 * NDATA;
        m_data[index + INDEX_ORDER] = 0.0;
        m_data[index + INDEX_MAX_CORR] = 0.0;
        m_data[index + INDEX_MAX_CORR_DT] = 0.0;

        let dCorrMax = -999999.0;
        let iCorrMax = -1;
        let iXShift = 0;
        let c = i0 * nbpl;
        if (i0 <= nCorr - 1 - iXShift || i0 >= m_nData - iXShift) {
            for (let k = 0, mCorrInd = i0; k < ih; k++, mCorrInd += m_nData, c += lw) {
                m_image.data[c] = 0;
                m_image.data[c + 1] = 0;
                m_image.data[c + 2] = 0;
            }
            return 0.0;
        }
        let dCorrScale = (m_maxCorr - m_minCorr) > 0.0 ? 1.0 / (m_maxCorr - m_minCorr) : 0.0;
        let i1Start = i0 - nCorr;
        let i1Min = nCorr - iXShift;
        let i1Max = m_nData - iXShift - 1;
        let ig = HIGHLIGHT;
        let kCorrMax = -1;
        
        let sOrder = 0.0;
        let sCorr = 0.0;
        let dMaxCorr = m_maxCorr - m_minCorr;
        if (dMaxCorr > 0.0001)
            dMaxCorr = 1.0 / dMaxCorr;
        let yCorr = 0.0;
        let dOrder = 0.0;
        for (let k = 0, iCorrInd = i0; k < ih; k++, iCorrInd += m_nData) {
            let ind = i1Start - k;
            if (ind >= i1Min && ind <= i1Max) {

                let corr = m_Corr[iCorrInd - 1];
                /*
                let indn = ind * NDATA;
                indn += NDATA;
                dOrder = m_data[indn + INDEX_HIGH] - m_data[indn + INDEX_OPEN];
                dOrder -= sg_main.m_Spread;
                if (dOrder > 0.0) {
                    dOrder = (m_data[indn + INDEX_OPEN] > 0.0) ? 1.0 / m_data[indn + INDEX_OPEN] : 0.0;
                    if (m_maxCorr > m_minCorr) {
                        let dCorr = (corr - m_minCorr) * dMaxCorr;
                        //dCorr *= dCorr;
                        dOrder = dOrder * dCorr * dCorr;
                        sOrder += dOrder;
                        sCorr += dCorr;
                        yCorr += k * dCorr;
                    }
                }
                else
                    dOrder = 0.0;
                */
                
                if (corr > dCorrMax) {
                    dCorrMax = corr;
                    kCorrMax = k;
                }
                
            }
        }
        if (kCorrMax < 0) {
            for (let k = 0, mCorrInd = i0; k < ih; k++, mCorrInd += m_nData, c += lw) {
                m_image.data[c] = 0;
                m_image.data[c + 1] = 0;
                m_image.data[c + 2] = 0;
            }
            return 0.0;
        }

        if (m_prevMaxCorrPos < 0 && kCorrMax >= 0) {
            m_prevMaxCorrPos = kCorrMax;
            m_prevMaxCorr = dCorrMax;
        }

        for (let k = 0, mCorrInd = i0; k < ih; k++, mCorrInd += m_nData, c += lw)
        {
            let i1 = i1Start - k;
            if (i1 >= i1Min && i1 <= i1Max)
            {
                let ir;
                let ib;
                let corrToDraw = m_Corr[mCorrInd - 1];
                let ccorr = (corrToDraw - m_minCorr) * dCorrScale;
                if (ccorr > 0.5)
                {
                    let rr = (ccorr - 0.5) * 510 + HIGHLIGHT;
                    ir = rr > 255.0 ? 255 : Math.floor(rr);
                    ib = 0;
                }
                else
                {
                    let bb = (0.5 - ccorr) * 510 + HIGHLIGHT;
                    ib = (bb > 255.0) ? 255 : Math.floor(bb);
                    ir = 0;
                }
                m_image.data[c] = ir;
                m_image.data[c+1] = ig;
                m_image.data[c+2] = ib;
            }
            else
            {
                m_image.data[c] = 0;
                m_image.data[c + 1] = 0;
                m_image.data[c + 2] = 0;
            }
        }
        dOrder = getOrder(dCorrMax, i0, nCorr, kCorrMax);
        m_data[index + INDEX_ORDER] = dOrder;
        m_data[index + INDEX_MAX_CORR] = dCorrMax;
        m_data[index + INDEX_MAX_CORR_DT] = kCorrMax;

        return dOrder;
    }
    histDataToTable = function (str) {
        if (localStorage != undefined) {
            localStorage.setItem("sg_hist_data", str);
        }

        let histTable = document.getElementById("historical_data");
        let allRows = str.split(/\r?\n|\r/);
        let table = "<table id='historical_data_table'>";
        let indData = 0;
        m_nData = 0;
        for (let singleRow = 0; singleRow < allRows.length; singleRow++) {
            if (singleRow === 0) {
                table += '<thead>';
                table += '<tr>';
            }
            else {
                table += '<tr>';
            }
            let rowCells = allRows[singleRow].split(',');
            if (rowCells.length >= 5) {
                for (let rowCell = 0; rowCell < rowCells.length; rowCell++) {
                    if (singleRow === 0) {
                        table += '<th>';
                        table += rowCells[rowCell];
                        table += '</th>';
                    }
                    else {
                        table += '<td>';
                        table += rowCells[rowCell];
                        table += '</td>';
                    }
                    if (singleRow != 0) {
                        if (rowCell == 0) {
                            m_DateTime[m_nData] = rowCells[rowCell];
                            m_nData++;
                        }
                        else
                            m_data[indData + rowCell - 1] = parseFloat(rowCells[rowCell]);
                    }
                }
                if (singleRow != 0) {
                    for (let fillRaw = rowCells.length - 1; fillRaw < NDATA; fillRaw++) {
                        m_data[indData + fillRaw] = 0.0;
                    }
                    indData += NDATA;
                }
            }

            if (singleRow === 0) {
                table += '</tr>';
                table += '</thead>';
                table += '<tbody>';
            }
            else {
                table += '</tr>';
            }
        }
        table += '</tbody>';
        table += '</table>';
        histTable.innerHTML = table;
    }

    loadData = function(myFile)
    {
        dt3 = new Date();

        m_Start = 0;
        m_Stop = 0;
        if (myFile)
        {
            let reader = new FileReader();
            reader.onloadend = function (e) { dataLoaded() };
            m_nData = 0;
            reader.addEventListener('load', function (e) { histDataToTable(e.target.result); });
            reader.readAsText(myFile);
        }
    }
    annotate = function()
    {
        if (m_Index0 >= m_Start && m_Index0 < m_Stop)
        {
            for (let i = 0; i < m_gDraw.length; i++)
            {
                if (m_gDraw[i])
                {
                    let texth = 12;
                    let texta = 4;
                    let v = m_data[m_Index0 * NDATA + i];
                    let dataScale = 1.0;
                    if (i == INDEX_ORDER) {
                        v = v * sg_main.m_OrderScale * m_data[m_Index0 * NDATA + INDEX_HIGH];
                    }
                    let x = getInterpolatedValue(m_Index0, m_imWidth, m_Start, m_Stop);
                    let y = m_imHeight - (v - m_gMin[i]) * m_imHeight / (m_gMax[i] - m_gMin[i]) - 1;
                    let txt = m_gName[i] + " " + parseFloat(v).toFixed(4);
                    if (y < texth / 2 + 2)
                        y = texth / 2 + 2;
                    else if (y > m_imHeight - texth / 2)
                        y = m_imHeight - texth / 2;
                    else if (x > m_imWidth / 2)
                        txt = txt + " -";
                    else
                        txt = "-" + txt;

                    sg_main.m_ctx.font = texth + "pt sans-serif";
                    let wh = sg_main.m_ctx.measureText(txt);
                    if (x > m_imWidth / 2) {
                        x = x - wh.width;
                    }
                    sg_main.m_ctx.fillStyle = "rgba(0, 0, 0, 0.4)";
                    sg_main.m_ctx.fillRect(x, y - texth / 2 - texta, wh.width, texth + texta * 2);
                    sg_main.m_ctx.fillStyle = m_gColors[i];
                    sg_main.m_ctx.fillText(txt, x, y + texth / 2);
                }
            }
        }
    }
    myDrawImage = function()
    {
        sg_main.m_ctx.lineCap = "butt";
        dt5 = new Date();
        dt6 = new Date();
        fillBitmap();
        dt7 = new Date();
        annotate();
        dt8 = new Date();

        if (onDrawEnd)
            onDrawEnd();

    }
    processLoadedData = function() {
        deleteCorrArrays();
        if (m_Stop == 0) {
            m_Start = 0;
            m_Stop = m_nData;
        }
        calculateArraysForCorrelation();
        myDrawImage();
    }
    dataLoaded = function () {
        dt4 = new Date();
        processLoadedData();
    }
    setup = function(x, y)
    {
        dt1 = new Date();
        m_imWidth = x - sg_main.m_cnv.offsetLeft*2;
        sg_main.m_cnv.width = m_imWidth;
        m_imHeight = y - sg_main.m_cnv.offsetTop;
        if (m_imHeight < 300)
            m_imHeight = 300;
        sg_main.m_cnv.height = m_imHeight;
        m_imageCanvas = sg_main.m_ctx.getImageData(0, 0, sg_main.m_cnv.width, sg_main.m_cnv.height);
        dt2 = new Date();
    }
    set1stOptimizeRanges = function()
    {
        sg_main.getOptimizeRanges();
    }
    run1stOptimization = function () {
        sg_main.m_OrderScale = m_ranges[0];
        sg_main.m_MaxLost = m_ranges[3];
        sg_main.m_AverageLength = m_ranges[6];
        sg_main.m_CorrLength = m_ranges[9];
        sg_main.m_CorrLength = m_ranges[10]; // go from max to min
        deleteCorrArrays();
        calculateArraysForCorrelation();
        setCorrelationMatrix();
    }
    set1stOptimizationValues = function () {
        m_optimum_PredScale = m_ranges[0];
        m_optimum_MaxLost = m_ranges[3];
        m_optimum_ave = m_ranges[6];
        m_optimum_ncorr = m_ranges[10]; // go from max to min
    }
    set2ndOptimization = function () {
        let optValues = [sg_main.m_OrderScale, sg_main.m_MaxLost, sg_main.m_AverageLength, sg_main.m_CorrLength];
        let iStep = [4, 4, 4, 2];
        let nSteps = 10;
        for (let i = 0; i < 4; i++)
        {
            let rmin = optValues[i] - m_ranges[i * 3 + 2] * nSteps;
            if (rmin < m_ranges[i * 3])
                rmin = m_ranges[i * 3];
            let rmax = optValues[i] + m_ranges[i * 3 + 2] * nSteps;
            if (rmax > m_ranges[i * 3 + 1])
                rmax = m_ranges[i * 3 + 1];
            m_ranges[i * 3] = rmin;
            m_ranges[i * 3 + 1] = rmax;
        }

        run1stOptimization();
    }
    set1stOptimization = function () {
        m_optimumProfit = -10000000000.0;

        set1stOptimizeRanges();
        set1stOptimizationValues();

        run1stOptimization();
    }
    setCurrentOptimizedValues = function () {
        sg_main.m_OrderScale = m_optimum_PredScale;
        sg_main.m_MaxLost = m_optimum_MaxLost;
        sg_main.m_AverageLength = m_optimum_ave;
        sg_main.m_CorrLength = m_optimum_ncorr;
    }
    stopOptimize = function () {
        setCurrentOptimizedValues();
        if (onOptimizeEnd)
            onOptimizeEnd();
    }
    optimizeParams = function() {
        if (m_optimizeTimeout)
            clearTimeout(m_optimizeTimeout);
        switch (m_optimize)
        {
            case eStart1stOptimizationLoop:
                set1stOptimization();
                calculateArraysForCorrelation();
                m_optimize++;
                break;
            case eRun1stOptimizationLoop:
                optomizeLoop(m_optimize);
                break;
            case eStop1stOptimizationLoop:
                setCurrentOptimizedValues();
                m_optimize++;
                break;
            case eStart2ndOptimizationLoop:
                set2ndOptimization();
                calculateArraysForCorrelation();
                m_optimize++;
                break;
            case eRun2ndOptimizationLoop:
                optomizeLoop(m_optimize);
                break;
            case eStop2ndOptimizationLoop:
                m_optimize++;
                break;
            case eEndOfOptimization:
                m_optimize = eNoOptimization;
                stopOptimize();
                break;
        }
        if (m_optimize != eNoOptimization)
        {
            sg_main.updateSlidersValues();
            m_optimizeTimeout = window.setTimeout(optimizeParams, 20);
        }
    }
    updateOptimums = function()
    {
        let s = getProfit();
        if (m_optimumProfit < s)
        {
            m_optimumProfit = s;
            m_optimum_PredScale = sg_main.m_OrderScale;
            m_optimum_MaxLost = sg_main.m_MaxLost;
            m_optimum_ave = sg_main.m_AverageLength;
            m_optimum_ncorr = sg_main.m_CorrLength;
            processLoadedData();
        }
    }
    optomizeLoop = function(iopt)
    {
        let iRanges = m_ranges;
        let iTick = new Date();
        let iStepCorr = (m_optimize == eRun2ndOptimizationLoop) ? 1 : 2;
        let iStep = (m_optimize == eRun2ndOptimizationLoop) ? 1 : 4;
        while(m_optimize == iopt)
        {
            updateOptimums();
            if (sg_main.m_OrderScale >= iRanges[1])
            {
                sg_main.m_OrderScale = iRanges[1];
                if (sg_main.m_OrderScale != iRanges[0]) // skip if iRanges[0] == iRanges[1]
                {
                    sg_main.m_OrderScale = iRanges[0];
                    calculateOrderDelta();
                }
                if (sg_main.m_MaxLost >= iRanges[4])
                {
                    sg_main.m_MaxLost = iRanges[4];
                    if (sg_main.m_MaxLost != iRanges[3]) // skip if iRanges[3] == iRanges[4]
                    {
                        sg_main.m_MaxLost = iRanges[3];
                        calculateOrderDelta();
                    }
                    if (sg_main.m_AverageLength >= iRanges[7])
                    {
                        sg_main.m_AverageLength = iRanges[7];
                        if (sg_main.m_AverageLength != iRanges[6])  // skip if iRanges[6] == iRanges[7]
                        {
                            sg_main.m_AverageLength = iRanges[6];
                            deleteCorrArrays();
                            setCorrelationMatrix();
                        }
                        if (sg_main.m_CorrLength <= iRanges[9]) // go from max to min
                        {
                            sg_main.m_CorrLength = iRanges[9]; // go from max to min
                            m_optimize++;
                            break;
                        }
                        else
                        {
                            sg_main.m_CorrLength -= iStepCorr * iRanges[11]; // go from max to min
                            if (sg_main.m_CorrLength <= iRanges[9]) // go from max to min
                                sg_main.m_CorrLength = iRanges[9]; // go from max to min
                            deleteCorrArrays();
                            setCorrelationMatrix();
                        }
                    }
                    else
                    {
                        sg_main.m_AverageLength += iStep * iRanges[8];
                        if (sg_main.m_AverageLength >= iRanges[7])
                            sg_main.m_AverageLength = iRanges[7];
                        deleteCorrArrays();
                        setCorrelationMatrix();
                    }
                }
                else
                {
                    sg_main.m_MaxLost += iStep * iRanges[5];
                    if (sg_main.m_MaxLost >= iRanges[4])
                        sg_main.m_MaxLost = iRanges[4];
                    calculateMaxLost();
                }
            }
            else
            {
                sg_main.m_OrderScale += iStep * iRanges[2];
                calculateOrderDelta();
            }
            let dt = new Date();
            let iDt = dt.getTime() - iTick.getTime();
            if (iDt > 200)
                break;
        }
    }
    startOrStopOptimization = function () {
        if (m_optimize == eNoOptimization) {
            m_optimize = eStart1stOptimizationLoop;
            optimizeParams();
            return true;
        }
        else {
            m_optimize = eStop2ndOptimizationLoop;
        }
        return false;
    }
    getPerfomance = function () {
        let strPerfomance =
            "Setup: " + (dt2.getTime() - dt1.getTime()).toString() + "[ms], " +
            "Data load: " + (dt4.getTime() - dt3.getTime()).toString() + "[ms], " +
            "Background setup: " + (dt6.getTime() - dt5.getTime()).toString() + "[ms], " +
            "Correlation map: " + (dt7.getTime() - dt6.getTime()).toString() + "[ms], " +
            "Draw to screen: " + (dt8.getTime() - dt7.getTime()).toString() + "[ms]<br>" +
            "Draw total: " + (dt8.getTime() - dt5.getTime()).toString() + "[ms] ";
        return strPerfomance;
    }
    return {
        set onOptimizeEnd(v) { onOptimizeEnd = v; },
        set onDrawEnd(v) { onDrawEnd = v; },
        get m_profit() { return m_profit; },
        get m_nData() { return m_nData; }, set m_nData(v) { m_nData = v; },
        get m_Start() { return m_Start; }, set m_Start(v) { m_Start = v; },
        get m_Stop() { return m_Stop; }, set m_Stop(v) { m_Stop = v; },
        get m_optimize() { return m_optimize; }, set m_optimize(v) { m_optimize = v; },
        get m_Index0() { return m_Index0; }, set m_Index0(v) { m_Index0 = v; },
        get m_Index1() { return m_Index1; }, set m_Index1(v) { m_Index1 = v; },
        setup: setup,
        myDrawImage : myDrawImage,
        processLoadedData: processLoadedData,
        zoom: zoom,
        orderString: orderString,
        loadData: loadData,
        getPerfomance: getPerfomance,
        startOrStopOptimization: startOrStopOptimization,
        m_DateTime: m_DateTime,
        m_data: m_data,
        m_ranges: m_ranges,
        histDataToTable: histDataToTable
    };
}
