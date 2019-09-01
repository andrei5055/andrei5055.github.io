// SGraph - Create bitmap with autocorrelation colormapData and optional graphics
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
    /* following "var" must be on top of all other attributes definitions*/
    let ind = 0;
    var INDEX_OPEN = ind; ind++;
    var INDEX_HIGH = ind; ind++;
    var INDEX_LOW = ind; ind++;
    var INDEX_CLOSE = ind; ind++;
    var INDEX_ADJ_CLOSE = ind; ind++;
    var INDEX_VOLUME = ind; ind++;
    var INDEX_PROFIT = ind; ind++;
    var INDEX_MAX_CORR = ind; ind++;
    var INDEX_MAX_CORR_DT = ind; ind++;
    var NDATA = ind;

    var INDEX_ORDER = NDATA; // Fake index: data are not in m_data array, but in array m_order

    ind = 0;
    var INDEX_OPEN_FP = ind; ind++;
    var INDEX_OPEN_SPREAD_FP = ind; ind++;
    var INDEX_HIGH_SPREAD_FP = ind; ind++;
    var INDEX_LOW_SPREAD_FP = ind; ind++;
    var INDEX_CLOSE_FP = ind; ind++;
    var INDEX_CLOSE_SPREAD_FP = ind; ind++;
    var INDEX_MAX_LOST_FP = ind; ind++;
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

    var HIGHLIGHT = 35;
    var MIN_LENGTH = 100;

    /* add new attributes below this line */
    var dt1 = new Date();
    var dt2 = new Date();
    var dt3 = new Date();
    var dt4 = new Date();
    var dt5 = new Date();
    var dt6 = new Date();
    var dt7 = new Date();
    var dt8 = new Date();
    var dt9 = new Date();
    var dt10 = new Date();
    var m_loop1Time = 0;
    var m_loop1Count = 0;

    var onOptimizeEnd = null;
    var onDrawEnd = null;
    var m_colorMap = null;
    var m_colorMapData = null;
    var m_imageCanvas = null;
    var m_imWidth = 0;
    var m_imHeight = 0;
    var m_corrWidth = 0;
    var m_corrHeight = 0;
    var m_Index0 = -1;
    var m_Index1 = -1;

    var m_Start = 0;
    var m_Stop = 0;
    var m_nData = 0;
    var m_ranges = new Array(0);
    var m_optimumProfit = 0.0;
    var m_optimize = eNoOptimization;

    var m_DateTime = new Array();

    var progress = 0;

    var m_xMouseDown = 0;
    var m_yMouseDown = 0;
    var m_order = new Array();
    var m_data = new Array();
    var m_dataFP = new Array();
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

    var m_optimizeTimeout = 0;
    var m_optimum_OrderScale = 0;
    var m_optimum_MaxLost = 0;
    var m_optimum_ncorr = 0;
    var m_profit = 0.0;
    var m_maxCorr = 0.0;
    var m_minCorr = 0.0;
    var m_gpu = null;
    var m_dataGPUSet = null;
    var m_dataGPU = null;
    var m_dataFPGPU = null;
    var m_dataFPGPUSet = null;
    var m_profitMatrixGPU = null;
    var m_getCorrelationMatrixGPU = null;
    var m_corrMatrixGPU = null;
    var m_spreadUsed = -999999.0;
    var m_corrLengthUsed = -999999.0;
    var m_bCalcOrder = true;
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
    drawCandlestick = function(ix1, iww, openp, high, low, closep)
    {
        let closeOpenMax = Math.max(closep, openp);
        let closeOpenMin = Math.min(closep, openp);
        let iy = Math.floor(high);
        let iy1 = Math.floor(closeOpenMax);
        let iww2 = Math.floor(iww / 2);
            
        let urRect = 0;
        let ugRect = 255;
        let ubRect = 0;
        if (closep < openp)
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

    getProfitOpt = function (N, iStepOrder, M, iStepLost) {
        let newOpt = false;
        let profiDt1 = new Date();
        let orderScale = new Array(N);
        for (let i = 0; i < N; i++) {
            orderScale[i] = m_ranges[0] + i * m_ranges[2] * iStepOrder;
            if (orderScale[i] > m_ranges[1]) {
                orderScale[i] = m_ranges[1];
            }
        }

        let maxLost = new Array(M);
        for (let i = 0; i < M; i++) {
            maxLost[i] = m_ranges[3] + i * m_ranges[5] * iStepLost;
            if (maxLost[i] > m_ranges[4]) {
                maxLost[i] = m_ranges[4];
            }
        }
        if (sg_main.m_bUseGPU && m_optimize != eNoOptimization) {
            if (!m_dataFPGPU) {
                m_dataFPGPU = m_dataFPGPUSet(m_dataFP);
            }
            let istart = Math.floor(m_Start);
            let istop = Math.floor(m_Stop);
            let ss = m_profitMatrixGPU(m_dataFPGPU, m_order, m_nData, orderScale, maxLost,
                sg_main.m_nDaysToKeepOpen, sg_main.m_Spread, istart, istop, NDATA_FP,
                INDEX_OPEN_FP, INDEX_CLOSE_FP, INDEX_MAX_LOST_FP,
                INDEX_OPEN_SPREAD_FP, INDEX_HIGH_SPREAD_FP, INDEX_LOW_SPREAD_FP, INDEX_CLOSE_SPREAD_FP);
            for (let i = 0; i < M; i++) {
                for (let j = 0; j < N; j++) {
                    if (m_optimumProfit < ss[i][j]) {
                        m_optimumProfit = ss[i][j];
                        m_optimum_OrderScale = orderScale[j];
                        m_optimum_MaxLost = maxLost[i];
                        newOpt = true;
                    }
                }
            }
        }
        else {
            for (let i = 0; i < M; i++) {
                for (let j = 0; j < N; j++) {
                    let s = getProfit(orderScale[j], maxLost[i]);
                    if (m_optimumProfit < s) {
                        m_optimumProfit = s;
                        m_optimum_OrderScale = orderScale[j];
                        m_optimum_MaxLost = maxLost[i];
                        newOpt = true;
                    }
                }
            }
        }
        let profiDt2 = new Date();
        m_loop1Time += profiDt2.getTime() - profiDt1.getTime();
        m_loop1Count += N * M;

        if (newOpt) {
            //console.log("Calculated max profit=" + m_optimumProfit);
            sg_main.m_OrderScale = m_optimum_OrderScale;
            sg_main.m_MaxLost = m_optimum_MaxLost;
            if (m_optimum_ncorr != sg_main.m_CorrLength) {
                m_optimum_ncorr = sg_main.m_CorrLength;
                if (!sg_main.m_bUseGPU)
                    processLoadedData();
            }
        }
    }
    setupGPU = function (nOpt) {
        let iStep = (nOpt == 2) ? 1 : 2;
        let N = Math.floor((m_ranges[1] - m_ranges[0]) / (iStep * m_ranges[2]) + 1.1);
        let M = Math.floor((m_ranges[4] - m_ranges[3]) / (iStep * m_ranges[5]) + 1.1);
        setGPUKernels(N, iStep, M, iStep, m_nData, sg_main.MATRIX_Y);
    }
    getCorrelationMatrix = function (nCorr) {
        if (m_corrLengthUsed != sg_main.m_CorrLength) {
            dt9 = new Date();
            m_corrLengthUsed = sg_main.m_CorrLength;
            m_minCorr = 1e+17;
            m_maxCorr = -1e+17;
            m_Corr.length = 0;
            if (sg_main.m_bUseGPU && m_optimize != eNoOptimization) {
                getCorrelationMatrixGPU(nCorr);
            }
            else {
                getCorrelationMatrixCPU(nCorr);
                //checkCorrelationMatrixCPU(nCorr);
            }
            if (Math.abs(m_minCorr - m_maxCorr) < 0.0001) {
                m_minCorr = -1.0;
                m_maxCorr = 1.0;
            }
            dt10 = new Date();
        }
    }
    getCorrelationMatrixGPU = function (nCorr) {
        if (!m_dataGPU) {
            m_dataGPU = m_dataGPUSet(m_data);
        }
        let gpuCorr = m_getCorrelationMatrixGPU(m_dataGPU, nCorr, NDATA, INDEX_HIGH);
        let iend = sg_main.MATRIX_Y;
        if (gpuCorr.length == iend) {
            m_Corr = new Array(m_nData * iend);
            let k = 0;
            for (let i = 0; i < iend; i++) {
                for (let j = 0; j < m_nData; j++, k++) {
                    corr = gpuCorr[i][j];
                    if (m_minCorr > corr)
                        m_minCorr = corr;
                    if (m_maxCorr < corr)
                        m_maxCorr = corr;
                    m_Corr[k] = corr;
                }
            }
        }

    }

    checkCorrelationMatrixCPU = function (ncorr) {
        let indCorr = 0;
        let ny = sg_main.MATRIX_Y;
        let ncorr2 = ncorr * 2;
        let ncorr2p1 = ncorr2 + 1;
        let ncorrh = ncorr / 2;
        let ncorrh1h = (ncorr - 1) / 2;
        let dncorr = 1.0 / ncorr;
        let d2ncorr = 1.0 / ncorrh / ncorrh1h;
        let xc = new Array(ncorr);
        let diffMax = -999999.0;
        for (let x = 0; x < m_nData; x++) {
            indCorr = x;
            if (indCorr < ncorr2p1) {
                for (let y = 0; y < ny; y++, indCorr += m_nData) {
                    m_Corr[indCorr] = 0.0;
                }
            }
            let indx = x * NDATA + INDEX_HIGH;
            let avex1 = 0.0;
            let avex2 = 0.0;
            let indxt = indx;
            for (let i = 0; i < ncorr; i++) {
                if (i >= ncorrh) {
                    avex1 += m_data[indxt];
                }
                else {
                    avex2 += m_data[indxt];
                }
                indxt -= NDATA;
            }
            let avx = (avex1 + avex2) * dncorr;
            let dtx = (avex2 - avex1) * d2ncorr;
            avx += dtx * ncorrh1h;
            indxt = indx;
            let avxt = avx;
            let sqx = 0.0;
            for (let j = 0; j < ncorr; j++) {
                let dx = (m_data[indxt] - avxt);
                sqx += dx * dx;
                avxt -= dtx;
                indxt -= NDATA;
            }
            indxt = indx;
            for (let k = 0; k < ncorr; k++) {
                xc[k] = (m_data[indxt] - avx);
                avx -= dtx;
                indxt -= NDATA;
            }
            for (let y = 0; y < ny; y++) {
                if (x - y < ncorr2p1) {
                    m_Corr[indCorr] = 0.0;
                }
                else {
                    let indy = indx - (ncorr + y) * NDATA;
                    let avey1 = 0.0;
                    let avey2 = 0.0;
                    let indyt = indy;
                    for (let i = 0; i < ncorr; i++) {
                        if (i >= ncorrh) {
                            avey1 += m_data[indyt];
                        }
                        else {
                            avey2 += m_data[indyt];
                        }
                        indyt -= NDATA;
                    }
                    let avy = (avey1 + avey2) * dncorr;
                    let dty = (avey2 - avey1) * d2ncorr;
                    let corr = 0.0;
                    avy += dty * ncorrh1h;
                    indxt = indx;
                    indyt = indy;
                    let avyt = avy;
                    let sqy = 0.0;
                    for (let j = 0; j < ncorr; j++) {
                        let dy = (m_data[indyt] - avyt);
                        sqy += dy * dy;
                        avyt -= dty;
                        indyt -= NDATA;
                    }
                    let d = sqx * sqy;
                    if (d > 0.000001)
                        d = 1.0 / Math.sqrt(d);
                    else
                        d = 0.0;
                    for (let k = 0; k < ncorr; k++) {
                        corr += xc[k] * (m_data[indy] - avy);
                        avy -= dty;
                        indy -= NDATA;
                    }
                    corr = corr * d;
                    let diff = Math.abs(m_Corr[indCorr] - corr);
                    if (diffMax < diff)
                        diffMax = diff;
                    //m_Corr[indCorr] = corr;

                    if (m_minCorr > corr)
                        m_minCorr = corr;
                    if (m_maxCorr < corr)
                        m_maxCorr = corr;
                }
                indCorr += m_nData;
            }
        }
        if (diffMax > 0.001)
            alert("Corr matrix difference=" + diffMax);
    }
    getCorrelationMatrixCPU = function (ncorr) {
        let indCorr = 0;
        let matrixh = sg_main.MATRIX_Y;
        let ncorr2p1 = ncorr * 2 + 1;
        let ncorrm1 = ncorr - 1;
        let ncorrh = ncorr / 2;
        let ncorrhp1 = ncorrh + 1;
        let ncorrm1h = (ncorr - 1) / 2;
        let dncorr = 1.0 / ncorr;
        let d2ncorr = 1.0 / ncorrh / ncorrm1h;
        let xc = new Array(ncorr);
        let e1 = new Array(m_nData);
        let e2 = new Array(m_nData);
        let e3 = new Array(m_nData);
        m_Corr = new Array(m_nData * matrixh);
        let e1sum = 0.0;
        let e2sum = 0.0;
        let e3sum = 0.0;
        let e4sum = 0.0;
        for (let i = 0, indv = INDEX_HIGH; i < m_nData; i++, indv += NDATA) {
            let v = m_data[indv];
            e1sum += v;
            e2sum += v * v;
            e1[i] = e1sum;
            e2[i] = e2sum;
            
            if (i < ncorr) {
                e3sum += v * i;
                e4sum += i * i;
            }
            else {
               e3sum = e3sum + (v * (ncorr - 1) - (e1[i - 1] - e1[i - ncorr]));
            }
            e3[i] = e3sum;
        }
        for (let indCorr1 = 0; indCorr1 < m_nData; indCorr1++) {
            indCorr = indCorr1;
            if (indCorr < ncorr2p1) {
                for (let j = 0; j < matrixh; j++, indCorr += m_nData) {
                    m_Corr[indCorr] = 0.0;
                }
            }
            let indx = indCorr1 * NDATA + INDEX_HIGH;
            let avex2 = e1[indCorr1] - e1[indCorr1 - ncorrh];
            let avex1 = e1[indCorr1 - ncorrh] - e1[indCorr1 - ncorr];
            let avx = (avex1 + avex2) * dncorr;
            let dtx = (avex2 - avex1) * d2ncorr;
            let a = dtx;
            let b = avx - dtx * ncorrm1h;
            let sqx = e2[indCorr1] - e2[indCorr1 - ncorr] -
                2.0 * (e3[indCorr1] * a + (e1[indCorr1] - e1[indCorr1 - ncorr]) * b) +
                e4sum * a * a + ncorr * ncorrm1 * a * b + b * b * ncorr;
            avx += dtx * ncorrm1h;
            indxt = indx;
            for (let k = 0; k < ncorr; k++) {
                xc[k] = (m_data[indxt] - avx);
                avx -= dtx;
                indxt -= NDATA;
            }
            for (let m = 0; m < matrixh; m++) {
                if (indCorr1 - m < ncorr2p1) {
                    m_Corr[indCorr] = 0.0;
                }
                else {
                    let indCorr0 = indCorr1 - ncorr - m;
                    let indy = indx - (ncorr + m) * NDATA;
                    let avey2 = e1[indCorr0] - e1[indCorr0 - ncorrh];
                    let avey1 = e1[indCorr0 - ncorrh] - e1[indCorr0 - ncorr];
                    let avy = (avey1 + avey2) * dncorr;
                    let dty = (avey2 - avey1) * d2ncorr;
                    let a = dty;
                    let b = avy - dty * ncorrm1h;
                    let sqy = e2[indCorr0] - e2[indCorr0 - ncorr] -
                        2.0 * (e3[indCorr0] * a + (e1[indCorr0] - e1[indCorr0 - ncorr]) * b) +
                        e4sum * a * a + ncorr*ncorrm1*a*b + b*b*ncorr;
                    let d = sqx * sqy;
                    if (d > 0.000001)
                        d = 1.0 / Math.sqrt(d);
                    else
                        d = 0.0;
                    avy += dty * ncorrm1h;
                    let corr = 0.0;
                    for (let k = 0; k < ncorr; k++) {
                        corr += xc[k] * (m_data[indy] - avy);
                        avy -= dty;
                        indy -= NDATA;
                    }
                    corr = corr * d;
                    m_Corr[indCorr] = corr;

                    if (m_minCorr > corr)
                        m_minCorr = corr;
                    if (m_maxCorr < corr)
                        m_maxCorr = corr;
                }
                indCorr += m_nData;
            }
        }
    }
    setGPUKernels = function (N, iStepOrder, M, iStepLost, nData, nY) {
        m_dataGPUSet = null;
        m_dataFPGPUSet = null;
        m_profitMatrixGPU = null;
        m_getCorrelationMatrixGPU = null;
        if (m_dataGPUSet)
            m_dataGPUSet.destroy();
        if (m_dataFPGPUSet)
            m_dataFPGPUSet.destroy();
        if (m_profitMatrixGPU)
            m_profitMatrixGPU.destroy();
        if (m_getCorrelationMatrixGPU)
            m_getCorrelationMatrixGPU.destroy();
        if (m_gpu)
            m_gpu.destroy();
        let optData = {
            setPipeline: true,
            output: [nData * NDATA]
        };
        let optDataFP = {
            setPipeline: true,
            output: [nData * NDATA_FP]
        };
        let optProfit = {
            output: [N, M]
        };
        let optCorr = {
            output: [nData, nY]
        };
        if (typeof GPU == "undefined") {
            alert("GPU not defined. Check that you included support for GPU.\nGPU usage disabled");
            sg_main.m_bUseGPU = false;
            return;
        }
        if (!GPU || !GPU.isGPUSupported || !GPU.isSinglePrecisionSupported) {
            alert("GPU not supported on your PC or by this Browser!\nGPU usage disabled");
            sg_main.m_bUseGPU = false;
            return;
        }
        if (!GPU.isSinglePrecisionSupported) {
            alert("GPU floating point operations not supported on your PC or by this Browser!\nGPU usage disabled");
            sg_main.m_bUseGPU = false;
            return;
        }
        m_gpu = new GPU({ mode: 'gpu' });
        //m_gpu = new GPU({ mode: 'dev' });
        if (!m_gpu) {
            alert("Can't create GPU modulle\nGPU usage disabled");
            sg_main.m_bUseGPU = false;
            return;
        }
        m_dataGPUSet = m_gpu.createKernel(function (m_data) {
            return m_data[this.thread.x];
        }, optData);

        m_dataFPGPUSet = m_gpu.createKernel(function (m_dataFP) {
            return m_dataFP[this.thread.x];
        }, optDataFP);

        m_getCorrelationMatrixGPU = m_gpu.createKernel(function (m_data, ncorr, NDATA, INDEX_HIGH) {
            let ncorr2 = ncorr * 2;
            if (this.thread.x - this.thread.y < ncorr2) {
                return 0.0;
            }
            else {
                let ncorrh = ncorr / 2;
                let ncorrm1h = (ncorr - 1) / 2;
                let dncorr = 1.0 / ncorr;
                let d2ncorr = 1.0 / ncorrh / ncorrm1h;
                let indx = this.thread.x * NDATA + INDEX_HIGH;
                let indy = indx - (ncorr + this.thread.y) * NDATA;
                let avex1 = 0.0;
                let avey1 = 0.0;
                let avex2 = 0.0;
                let avey2 = 0.0;
                let indxt = indx;
                let indyt = indy;
                for (let i = 0; i < ncorr; i++) {
                    if (i >= ncorrh) {
                        avex1 += m_data[indxt];
                        avey1 += m_data[indyt];
                    }
                    else {
                        avex2 += m_data[indxt];
                        avey2 += m_data[indyt];
                    }
                    indxt -= NDATA;
                    indyt -= NDATA;
                }
                let avx = (avex1 + avex2) * dncorr;
                let avy = (avey1 + avey2) * dncorr;
                let dtx = (avex2 - avex1) * d2ncorr;
                let dty = (avey2 - avey1) * d2ncorr;
                avx += dtx * ncorrm1h;
                avy += dty * ncorrm1h;
                indxt = indx;
                indyt = indy;
                let avxt = avx;
                let avyt = avy;
                let sqx = 0.0;
                let sqy = 0.0;
                for (let j = 0; j < ncorr; j++) {
                    let dx = (m_data[indxt] - avxt);
                    sqx += dx * dx;
                    let dy = (m_data[indyt] - avyt);
                    sqy += dy * dy;
                    avxt -= dtx;
                    avyt -= dty;
                    indxt -= NDATA;
                    indyt -= NDATA;
                }
                let d = sqx * sqy;
                if (d > 0.000001)
                    d = 1.0 / Math.sqrt(d);
                else
                    d = 0.0;
                let s = 0;
                for (let k = 0; k < ncorr; k++) {
                    s += (m_data[indx] - avx) * (m_data[indy] - avy);
                    avx -= dtx;
                    avy -= dty;
                    indx -= NDATA;
                    indy -= NDATA;
                }
                s = s * d;
                return s;
            }
        }, optCorr);
        m_profitMatrixGPU = m_gpu.createKernel(function (m_dataFP, m_order, nData, orderScale, maxLost,
                nDaysToKeepOpen, spread, iStartP, iStopP, NDATA_FP,
                INDEX_OPEN_FP, INDEX_CLOSE_FP, INDEX_MAX_LOST_FP,
                INDEX_OPEN_SPREAD_FP, INDEX_HIGH_SPREAD_FP, INDEX_LOW_SPREAD_FP, INDEX_CLOSE_SPREAD_FP) {
            let s = 0.0;
            let buyPrice = 0.0;
            let bOrder = false;
            let sellLimit = buyPrice;
            let prof = 0.0;
            let closep = 0.0;
            let todayProfit = 0.0;
            let nDaysOpenPos = 0;
            let iStart = iStartP;
            let iStop = iStopP;
            if (iStart < 2)
                iStart = 2;
            if (iStop < nData)
                iStop++;
            let vMaxLost = maxLost[this.thread.y];
            let vOrderScale = orderScale[this.thread.x];
            let indData = iStart * NDATA_FP;
            for (let i = 0; i < nData; i++) {
                if (i >= iStop) {
                    break;
                }
                else {
                    if (i >= iStart) {
                        let openp = m_dataFP[indData + INDEX_OPEN_FP];
                        let highSpread = m_dataFP[indData + INDEX_HIGH_SPREAD_FP];
                        closep = m_dataFP[indData + INDEX_CLOSE_FP];
                        let dMaxLost = m_dataFP[indData + INDEX_MAX_LOST_FP] * vMaxLost;
                        let delta = m_order[i] * vOrderScale;

                        todayProfit = 0.0;

                        if (!bOrder) {
                            if (delta > spread) {
                                buyPrice = openp + spread;
                                sellLimit = openp + delta;
                                bOrder = true;
                                nDaysOpenPos = 0;
                            }
                        }

                        if (bOrder) {
                            if (dMaxLost != 0.0) {
                                let openSpread = m_dataFP[indData + INDEX_OPEN_SPREAD_FP];
                                if (buyPrice - openSpread > dMaxLost) {
                                    todayProfit = openSpread - buyPrice;
                                    bOrder = false;
                                }
                                else {
                                    if ((sellLimit <= highSpread) && (openp > closep)) {
                                        // Assume (if open > close) sell limit triggered before lost limit
                                        bOrder = true; // Nop
                                    }
                                    else if (buyPrice - m_dataFP[indData + INDEX_LOW_SPREAD_FP] > dMaxLost) {
                                        todayProfit = -dMaxLost;
                                        bOrder = false;
                                    }
                                }
                            }
                            if (bOrder) {
                                if (sellLimit <= highSpread) {
                                    todayProfit = sellLimit - buyPrice;
                                    bOrder = false;
                                }
                                else {
                                    if (nDaysOpenPos < nDaysToKeepOpen) {
                                        nDaysOpenPos++;
                                    }
                                    else {
                                        prof = m_dataFP[indData + INDEX_CLOSE_SPREAD_FP] - buyPrice;
                                        // close position at the end of the day
                                        todayProfit = prof;
                                        bOrder = false;
                                    }
                                }
                            }
                        }
                        s += todayProfit;
                        indData += NDATA_FP;
                    }
                }
            }
            if ((iStop == nData) && bOrder) {
                todayProfit = closep - spread - buyPrice;
                s += todayProfit;
            }
           // return s;
            return indData + s / 100.0;
        }, optProfit);

        if (!m_dataGPUSet || !m_dataFPGPUSet || !m_getCorrelationMatrixGPU || !m_profitMatrixGPU) {
            alert("Error in one of the GPU function\nGPU usage disabled");
            m_dataGPUSet = null;
            m_dataFPGPUSet = null;
            m_profitMatrixGPU = null;
            m_getCorrelationMatrixGPU = null;
            sg_main.m_bUseGPU = false;
            return;
        }
        m_profitMatrixGPU.setLoopMaxIterations(m_nData);
    }
    getProfit = function(vOrderScale, vMaxLost)
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
        let closep = 0.0;
        let todayProfit = 0.0;
        let nDaysOpenPos = 0;
        let nDaysToKeepOpen = sg_main.m_nDaysToKeepOpen;
        let iStart = Math.floor(m_Start);
        let iStop = Math.floor(m_Stop);
        if (iStart < 2)
            iStart = 2;
        if (iStop < m_nData)
            iStop++;
        for (let i = iStart, indData = iStart * NDATA_FP, outData = iStart * NDATA; i < iStop; i++, indData += NDATA_FP, outData += NDATA)
        {
            let openp = m_dataFP[indData + INDEX_OPEN_FP];
            let highSpread = m_dataFP[indData + INDEX_HIGH_SPREAD_FP];
            closep = m_dataFP[indData + INDEX_CLOSE_FP];
            let dMaxLost = m_dataFP[indData + INDEX_MAX_LOST_FP] * vMaxLost;
            let delta = m_order[i] * vOrderScale;
  
            todayProfit = 0.0;
        
            if (!bOrder)
            {
                if (delta > spread)
                {
                    buyPrice = openp + spread;
                    sellLimit = openp + delta;
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
                        if (sellLimit <= highSpread && openp > closep)
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
            m_data[outData + INDEX_PROFIT] = s;
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
            todayProfit = closep - spread - buyPrice;
            s += todayProfit;
        }
        m_dMin = dMin;
        m_dMax = dMax;
        m_dMinPr = dMinPr;
        m_dMaxPr = dMaxPr;
        m_profit = s;

        //console.log("Calculated dis profit=" + s);
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
                let openp  = m_data[j + INDEX_OPEN];
                let high  = m_data[j + INDEX_HIGH];
                let low   = m_data[j + INDEX_LOW];
                let closep = m_data[j + INDEX_CLOSE];
                openp = ((openp - dMin) * dScale);
                high = ((high - dMin) * dScale);
                low = ((low - dMin) * dScale);
                closep = ((closep - dMin) * dScale);
                drawCandlestick(ix1, iww, openp, high, low, closep);
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
                v = m_order[i] * dataScale;
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
    scaleColormapToImage = function (colormapData, width, height, w, h, dStart, dStop) {
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
                    imageCanvas_data[cOut] = (colormapData[cI1] * (1.0 - yd) + colormapData[cI2] * yd);
                    imageCanvas_data[cOut + 1] = (colormapData[cI1 + 1] * (1.0 - yd) + colormapData[cI2 + 1] * yd);
                    imageCanvas_data[cOut + 2] = (colormapData[cI1 + 2] * (1.0 - yd) + colormapData[cI2 + 2] * yd);
                    imageCanvas_data[cOut + 3] = 255;
                }
                else {
                    imageCanvas_data[cOut] = ((colormapData[cI1] * (1.0 - xd) + colormapData[cI1 + imgBytesPerPixel] * xd) * (1 - yd) + (colormapData[cI2] * (1.0 - xd) + colormapData[cI2 + imgBytesPerPixel] * xd) * yd);
                    imageCanvas_data[cOut + 1] = ((colormapData[cI1 + 1] * (1.0 - xd) + colormapData[cI1 + (imgBytesPerPixel + 1)] * xd) * (1 - yd) + (colormapData[cI2 + 1] * (1.0 - xd) + colormapData[cI2 + (imgBytesPerPixel + 1)] * xd) * yd);
                    imageCanvas_data[cOut + 2] = ((colormapData[cI1 + 2] * (1.0 - xd) + colormapData[cI1 + (imgBytesPerPixel + 2)] * xd) * (1 - yd) + (colormapData[cI2 + 2] * (1.0 - xd) + colormapData[cI2 + (imgBytesPerPixel + 2)] * xd) * yd);
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
        setCorrelationMatrix();
        sg_main.m_ctx.clearRect(0, 0, m_imWidth, m_imHeight);
        scaleColormapToImage(m_colorMapData, iw, ih, m_colorMap.width, m_colorMap.height, m_Start, m_Stop);

        sg_main.m_ctx.putImageData(m_imageCanvas, 0, 0);

        drawGrpah("Open", ind, r, g, b, dMinA, dMaxA, sg_main.m_bShowOpen || sg_main.m_bShowCandleStick);
   
        highlightCorrIntervals(nCorr);
        
        let vOrderScale = sg_main.m_OrderScale;
        let vMaxLost = sg_main.m_MaxLost;

        getProfit(vOrderScale, vMaxLost);

        let dd = m_dMax - m_dMin;
        let dMin = m_dMin - dd * 0.5;
        let dMax = m_dMax + dd * 0.5;
        drawGrpah("Volume",      INDEX_VOLUME,      200, 200, 200, volumeMinMax.dMin,  volumeMinMax.dMax,    sg_main.m_bShowVolume);
        drawGrpah("Order",       INDEX_ORDER,       127, 127, 255, dMin,   dMax,     sg_main.m_bShowOrder);
        drawGrpah("Profit",      INDEX_PROFIT,      255, 255, 255, m_dMinPr, m_dMaxPr,   sg_main.m_bShowProfit);
        drawGrpah("Correlation", INDEX_MAX_CORR,    255, 127, 127,   -1.0,   1.0,      sg_main.m_bShowCorr);
        drawGrpah("Maximum Correlation", INDEX_MAX_CORR_DT, 255, 0, 255, 0.0, sg_main.MATRIX_Y, sg_main.m_bShowCorrPos);

        return 1;
    }
    calculateMaxLost = function () {
        for (let i = 0, indData = 0, outdData = 0; i < m_nData; i++, indData += NDATA, outdData += NDATA_FP) {
            m_dataFP[outdData + INDEX_MAX_LOST_FP] = m_data[indData + INDEX_OPEN] * 0.01;
        }
    }
    getAdjustedBySpreadOHLC = function()
    {
        let spread = sg_main.m_Spread;
        for (let i = 0, indData = 0, outdData = 0; i < m_nData; i++, indData += NDATA, outdData += NDATA_FP)
        {
            let openp = m_data[indData + INDEX_OPEN];
            let high = m_data[indData + INDEX_HIGH];
            let low = m_data[indData + INDEX_LOW];
            let closep = m_data[indData + INDEX_CLOSE];
            let openSpread  = openp  - spread;
            let highSpread  = high  - spread;
            let lowSpread   = low   - spread;
            let closeSpread = closep - spread;
            m_dataFP[outdData + INDEX_OPEN_FP]         = openp;
            m_dataFP[outdData + INDEX_CLOSE_FP]        = closep;
            m_dataFP[outdData + INDEX_OPEN_SPREAD_FP]  = openSpread;
            m_dataFP[outdData + INDEX_HIGH_SPREAD_FP]  = highSpread;
            m_dataFP[outdData + INDEX_LOW_SPREAD_FP]   = lowSpread;
            m_dataFP[outdData + INDEX_CLOSE_SPREAD_FP] = closeSpread;
        }
    }
    getArraysForProfitCalc = function()
    {
        if (Math.abs(sg_main.m_Spread - m_spreadUsed) > 0.00001) {
            m_spreadUsed = sg_main.m_Spread;
        getAdjustedBySpreadOHLC();
        calculateMaxLost();
            m_dataFPGPU = null;
        }
    }
    deleteCorrArrays = function()
    {
        m_bCalcOrder = true;
    }
    allocateImage = function(nx, ny)
    {
        if (m_colorMap == null || m_corrWidth != nx || m_corrHeight != ny) {
            m_corrWidth = nx;
            m_corrHeight = ny;
            return sg_main.m_ctx.getImageData(0, 0, nx, ny);
        }
        return m_colorMap;
    }
    calcColormapAndOrders = function (nCorr) {
        for (let i = 1; i < m_nData; i++) {
            calcOneOrderAndColormapColumn(m_nData - i, nCorr);
        }
        // Set first value
        m_order[0] = m_order[1];
        m_data[INDEX_MAX_CORR] = m_data[NDATA + INDEX_MAX_CORR];
        m_data[INDEX_MAX_CORR_DT] = m_data[NDATA + INDEX_MAX_CORR_DT];
    }
    setCorrelationMatrix = function()
    { 
        if (m_bCalcOrder || m_corrLengthUsed != sg_main.m_CorrLength)
        {
            let nCorr = sg_main.m_CorrLength * 2;
            if (nCorr < 2)
                nCorr = 2;

        m_colorMap = allocateImage(m_nData, sg_main.MATRIX_Y);
        m_colorMapData = m_colorMap.data;

            getArraysForProfitCalc();
            getCorrelationMatrix(nCorr);

            calcColormapAndOrders(nCorr);
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
    function setTo0colormap(i0) {
        // Set column of colormapData to 0
        let nbpl = 4;
        let ih = m_colorMap.height;
        let lw = m_colorMap.width * nbpl;
        let c = i0 * nbpl;
        for (let k = 0; k < ih; k++, c += lw) {
            m_colorMapData[c] = 0;
            m_colorMapData[c + 1] = 0;
            m_colorMapData[c + 2] = 0;
        }
    }
    function calcColormapColumn(i0, nCorr, i1Min, i1Max, dCorrScale) {
        // Set column of colormapData
        let nbpl = 4;
        let ih = m_colorMap.height;
        let lw = m_colorMap.width*nbpl;
        let ig = HIGHLIGHT;
        let c = i0 * nbpl;
        let i1Start = i0 - nCorr;
        for (let k = 0, mCorrInd = i0; k < ih; k++, mCorrInd += m_nData, c += lw) {
            let i1 = i1Start - k;
            if (i1 >= i1Min && i1 <= i1Max) {
                let ir;
                let ib;
                let corrToDraw = m_Corr[mCorrInd - 1];
                let ccorr = (corrToDraw - m_minCorr) * dCorrScale;
                if (ccorr > 0.5) {
                    let rr = (ccorr - 0.5) * 510 + HIGHLIGHT;
                    ir = rr > 255.0 ? 255 : Math.floor(rr);
                    ib = 0;
                }
                else {
                    let bb = (0.5 - ccorr) * 510 + HIGHLIGHT;
                    ib = (bb > 255.0) ? 255 : Math.floor(bb);
                    ir = 0;
                }
                m_colorMapData[c] = ir;
                m_colorMapData[c + 1] = ig;
                m_colorMapData[c + 2] = ib;
            }
            else {
                m_colorMapData[c] = 0;
                m_colorMapData[c + 1] = 0;
                m_colorMapData[c + 2] = 0;
            }
        }
    }
    calcOneOrderAndColormapColumn = function (i0, nCorr)
    {
        let index = i0 * NDATA;
        m_order[i0] = 0.0;
        m_data[index + INDEX_MAX_CORR] = 0.0;
        m_data[index + INDEX_MAX_CORR_DT] = 0.0;

        let iXShift = 0;
        if (i0 <= nCorr - 1 - iXShift || i0 >= m_nData - iXShift) {
            setTo0colormap(i0, nCorr);
            return 0.0;
        }
        let dCorrMax = -999999.0;
        let iCorrMax = -1;
        let dCorrScale = (m_maxCorr - m_minCorr) > 0.0 ? 1.0 / (m_maxCorr - m_minCorr) : 0.0;
        let i1Start = i0 - nCorr;
        let i1Min = nCorr - iXShift;
        let i1Max = m_nData - iXShift - 1;
        let kCorrMax = -1;
        
        let sOrder = 0.0;
        let sCorr = 0.0;
        let yCorr = 0.0;
        let dOrder = 0.0;
        let ih = m_colorMap.height;
        for (let k = 0, iCorrInd = i0; k < ih; k++, iCorrInd += m_nData) {
            let ind = i1Start - k;
            if (ind >= i1Min && ind <= i1Max) {
                let corr = m_Corr[iCorrInd - 1];
                if (corr > dCorrMax) {
                    dCorrMax = corr;
                    kCorrMax = k;
                }
            }
        }
        if (kCorrMax < 0) {
            setTo0colormap(i0);
            return 0.0;
        }
        dCorrMax = (dCorrMax - m_minCorr) * dCorrScale;
        calcColormapColumn(i0, nCorr, i1Min, i1Max, dCorrScale);
        dOrder = getOrder(dCorrMax, i0, nCorr, kCorrMax);
        m_order[i0] = dOrder * m_data[index + INDEX_OPEN];
        m_data[index + INDEX_MAX_CORR] = dCorrMax;
        m_data[index + INDEX_MAX_CORR_DT] = kCorrMax;
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
                        v = v * sg_main.m_OrderScale;
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
        myDrawImage();
    }
    dataLoaded = function () {
        dt4 = new Date();
        m_bCalcOrder = true;
        m_spreadUsed = -999999.0; // Recalculate internal arrays by getArraysForProfitCalc
        m_corrLengthUsed = -999999.0; // Recalculate correlation matrix by getCorrelationMatrix
        m_dataGPU = null; // resend loaded data to GPU (if 'Use GPU' selected)
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
    setupOptimization = function () {
        sg_main.m_OrderScale = m_ranges[0];
        sg_main.m_MaxLost = m_ranges[3];
        sg_main.m_CorrLength = m_ranges[7]; // go from max to min
        m_bCalcOrder = true;
        setCorrelationMatrix();
    }
    set1stOptimizeRanges = function()
    {
        sg_main.getOptimizeRanges();
    }
    set1stOptimizationValues = function () {
        m_optimum_OrderScale = m_ranges[0];
        m_optimum_MaxLost = m_ranges[3];
        m_optimum_ncorr = m_ranges[7]; // go from max to min
    }
    set2ndOptimization = function () {
        let optValues = [sg_main.m_OrderScale, sg_main.m_MaxLost, sg_main.m_CorrLength];
        //let nSteps = 10;
        let nSteps = 5;
        for (let i = 0; i < optValues.length; i++)
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

        if (sg_main.m_bUseGPU)
            setupGPU(2);

        setupOptimization();
    }
    set1stOptimization = function () {
        m_loop1Time = 0;
        m_loop1Count = 0;
        m_optimumProfit = -10000000000.0;

        set1stOptimizeRanges();
        set1stOptimizationValues();

        if (sg_main.m_bUseGPU)
            setupGPU(1);

        setupOptimization();
    }
    setCurrentOptimizedValues = function () {
        sg_main.m_OrderScale = m_optimum_OrderScale;
        sg_main.m_MaxLost = m_optimum_MaxLost;
        sg_main.m_CorrLength = m_optimum_ncorr;
    }
    stopOptimize = function () {
        setCurrentOptimizedValues();
        if (onOptimizeEnd)
            onOptimizeEnd();
    }
    optimizeParams = function() {
        if (m_optimizeTimeout) {
            clearTimeout(m_optimizeTimeout);
            m_optimizeTimeout = 0;
        }
        switch (m_optimize)
        {
            case eStart1stOptimizationLoop:
                set1stOptimization();
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
            let tm = 20;
            if (false || !!document.documentMode)
                tm = 100; // for IE delay need to be atleast 70 ms for sliders update 
            m_optimizeTimeout = window.setTimeout(optimizeParams, tm);
        }
    }
    updateOptimums = function (iStep) {
        let N = Math.floor((m_ranges[1] - m_ranges[0]) / (iStep * m_ranges[2]) + 1.1);
        let M = Math.floor((m_ranges[4] - m_ranges[3]) / (iStep * m_ranges[5]) + 1.1);
        getProfitOpt(N, iStep, M, iStep);
    }
    optomizeLoop = function (iopt) {
        let iRanges = m_ranges;
        let iTick = new Date();
        let iStepCorr = (m_optimize == eRun2ndOptimizationLoop) ? 1 : 2;
        let iStep = (m_optimize == eRun2ndOptimizationLoop) ? 1 : 2;
        while (m_optimize == iopt) {
            updateOptimums(iStep);
            if (sg_main.m_CorrLength == 2) // go from max to min
            {
                m_optimize++;
                break;
            }
            if (sg_main.m_CorrLength <= iRanges[6]) // go from max to min
            {
                sg_main.m_CorrLength = iRanges[6]; // go from max to min
                m_optimize++;
                break;
            }
            else {
                sg_main.m_CorrLength -= iStepCorr * iRanges[8]; // go from max to min
                if (sg_main.m_CorrLength <= iRanges[6]) // go from max to min
                    sg_main.m_CorrLength = iRanges[6]; // go from max to min
                if (sg_main.m_CorrLength < 2)
                    sg_main.m_CorrLength = 2;
                m_bCalcOrder = true;
                setCorrelationMatrix();
            }
            let dt = new Date();
            let iDt = dt.getTime() - iTick.getTime();
            if (iDt > 200)
                break;
        }
    }
    stopOptimization = function () {
        if (m_optimize != 0) {
            m_optimize = eStop2ndOptimizationLoop;
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
            "total 'OrderScale x MaxLost' loop=" + parseFloat(m_loop1Time/1000).toFixed(3) + "[sec], Average=" + parseFloat(m_loop1Time / (m_loop1Count+0.0001)).toFixed(3) + "[ms], " +
            "Colormap: " + (dt7.getTime() - dt6.getTime()).toString() + "[ms], " +
            "Corr Matrix: " + (dt10.getTime() - dt9.getTime()).toString() + "[ms], " +
            "Draw to screen: " + (dt8.getTime() - dt7.getTime()).toString() + "[ms]<br>" +
            "Draw total: " + (dt8.getTime() - dt5.getTime()).toString() + "[ms] ";
        
        //m_loop1Time = 0;
        //m_loop1Count = 0;
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
        stopOptimization: stopOptimization,
        startOrStopOptimization: startOrStopOptimization,
        m_DateTime: m_DateTime,
        m_data: m_data,
        m_ranges: m_ranges,
        histDataToTable: histDataToTable
    };
}
