#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <deque>
#include <vector>
#include <ctime>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <map>

const std::time_t WINDOW_DURATION_SECONDS = 60;
const std::time_t HIGH_LOW_WINDOW_DURATION = 120;
const std::time_t SECOND_WINDOW_DURATION = 45;
const std::time_t SEC_PER_PERIOD_MACD = 50;
const std::time_t SHORT_TERM_DURATION = SEC_PER_PERIOD_MACD * 12;  // 12 periods * 50 seconds
const std::time_t LONG_TERM_DURATION = SEC_PER_PERIOD_MACD * 26;  // 26 periods * 50 seconds
const std::time_t SIGNAL_LINE_DURATION = SEC_PER_PERIOD_MACD * 9; // 9 periods * 50 seconds
const int HISTORICAL_VOLATILITY_WINDOW = 30;
const double ANNUALIZATION_FACTOR = std::sqrt(252);
const int SMI_LOOKBACK_PERIOD = 600;  
const int SMI_SMOOTHING_PERIOD = 180; 
const int SMI_SIGNAL_PERIOD = 180; 
const int SMI_INTERVAL_SECONDS = 1;
const int RVI_PERIOD = 10; 
const int RVI_SIGNAL_PERIOD = 4; 
const int RSI_PERIOD = 14; 
const int RSI_SIGNAL_PERIOD = 9; 

struct TradeData {
    std::time_t timestamp;
    double volume;
    double totVol;
};

struct PriceData {
    std::time_t timestamp;
    double price;
};

struct ResampledPriceData {
    double high;
    double low;
    double close;
};

struct RVIData {
    double stdDevHigh = 0;
    double stdDevLow = 0;
    double stdDevClose = 0;
};

double totalBidVolume = 0;
double totalAskVolume = 0;
double windowHigh = 0;
double windowLow = std::numeric_limits<double>::max();
double periodVolume = 0;
double periodVolume2 = 0;

std::deque<TradeData> secondWindowDeque;
std::deque<double> macdLineDeque;  

double calculateSpread(double ask, double bid) {
    return ask - bid;
}

double calculateBookPressure(int bidSize, int askSize) {
    if (askSize == 0) {
        return 0;
    }
    return static_cast<double>(bidSize) / askSize;
}

double calculateAccumulationDistribution(double close, double low, double high, double periodVolume) {
    if (high == low) {
        return 0;
    }
    return (((close - low) - (high - close)) / (high - low)) * periodVolume;
}

double calculateHistoricalVolatility(const std::deque<PriceData>& prices) {
    if(prices.size() < 2) return 0;

    std::vector<double> logReturns;
    for(size_t i = 1; i < prices.size(); ++i) {
        double logReturn = std::log(prices[i].price / prices[i-1].price);
        logReturns.push_back(logReturn);
    }

    double sum = std::accumulate(logReturns.begin(), logReturns.end(), 0.0);
    double mean = sum / logReturns.size();

    double variance = 0.0;
    for(double lr : logReturns) {
        variance += (lr - mean) * (lr - mean);
    }
    variance /= (logReturns.size() - 1);

    double stdev = std::sqrt(variance);
    return stdev * ANNUALIZATION_FACTOR;
}

std::tm parseDateTime(const std::string& dateTime) {
    std::tm tm = {};
    std::istringstream ss(dateTime);
    ss >> std::get_time(&tm, "%m/%d/%y %H:%M:%S");
    return tm;
}

void updateHighLowDeque(std::deque<TradeData>& deque, std::time_t newTime, double& high, double& low) {
    while (!deque.empty() && newTime - deque.front().timestamp >= HIGH_LOW_WINDOW_DURATION) {
        deque.pop_front();
    }

    high = 0;
    low = std::numeric_limits<double>::max();
    for (const auto& trade : deque) {
        high = std::max(high, trade.volume);
        low = std::min(low, trade.volume);
    }
}


void updateDeque(std::deque<TradeData>& deque, std::time_t newTime, double& totalVolume) {
    while (!deque.empty() && newTime - deque.front().timestamp >= WINDOW_DURATION_SECONDS) {
        totalVolume -= deque.front().volume;
        deque.pop_front();
    }
}

void updateSecondWindow(std::deque<TradeData>& deque, std::time_t newTime) {
    while (!deque.empty() && newTime - deque.front().timestamp >= SECOND_WINDOW_DURATION) {
        deque.pop_front();
    }
}

void addDataToSecondWindow(std::time_t timestamp, double totVol) {
    secondWindowDeque.push_back({timestamp, 0, totVol});
    periodVolume2 += totVol;
}

double calculateTradingBalance() {
    if(totalBidVolume + totalAskVolume == 0) return 0;
    return totalBidVolume / (totalBidVolume + totalAskVolume);
}

double calculateEMA(const std::deque<PriceData>& prices, std::time_t duration) {
    if (prices.empty()) return 0;

    double multiplier = 2.0 / ((duration / SEC_PER_PERIOD_MACD) + 1); 
    double ema = prices.front().price;  // 

    for (const auto& priceData : prices) {
        if (priceData.timestamp - prices.front().timestamp <= duration) {
            ema = (priceData.price - ema) * multiplier + ema;
        } else {
            break;
        }
    }
    return ema;
}

double calculateEMAForMACD(const std::deque<double>& values, std::time_t duration) {
    if (values.empty()) return 0;

    double multiplier = 2.0 / ((duration / SEC_PER_PERIOD_MACD) + 1); 
    double ema = values.front();  

    for (size_t i = 1; i < values.size(); i++) {
        if (i * SEC_PER_PERIOD_MACD <= duration) {  
            ema = (values[i] - ema) * multiplier + ema;
        } else {
            break;
        }
    }
    return ema;
}

double getHighestHigh(const std::deque<PriceData>& prices) {
    double highestHigh = std::numeric_limits<double>::lowest();
    for (const auto& price : prices) {
        if (price.price > highestHigh) highestHigh = price.price;
    }
    return highestHigh;
}

double getLowestLow(const std::deque<PriceData>& prices) {
    double lowestLow = std::numeric_limits<double>::max();
    for (const auto& price : prices) {
        if (price.price < lowestLow) lowestLow = price.price;
    }
    return lowestLow;
}

double calculateEMA(const std::deque<double>& values, int periods) {
    if (values.empty()) return 0;
    double multiplier = 2.0 / (periods + 1);
    double ema = values.front();  // Starting with the first value

    for (size_t i = 1; i < values.size(); ++i) {
        ema = (values[i] - ema) * multiplier + ema;
    }
    return ema;
}

std::map<std::time_t, ResampledPriceData> resampleData(const std::deque<PriceData>& rawPrices) {
    std::map<std::time_t, ResampledPriceData> resampledData;
    for (const auto& price : rawPrices) {
        std::time_t bucket = price.timestamp - (price.timestamp % SMI_INTERVAL_SECONDS);
        auto& data = resampledData[bucket];
        if (data.high < price.price || data.high == 0) data.high = price.price;
        if (data.low > price.price || data.low == 0) data.low = price.price;
        data.close = price.price; 
    }
    return resampledData;
}

double standardDeviation(const std::deque<double>& data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0, std::plus<>(), [mean](double a, double b) { return (a - mean) * (b - mean); });
    return std::sqrt(sq_sum / data.size());
}

RVIData calculateRVIData(const std::deque<ResampledPriceData>& resampledData) {
    std::deque<double> highs, lows, closes;
    for (const auto& data : resampledData) {
        highs.push_back(data.high);
        lows.push_back(data.low);
        closes.push_back(data.close);
    }

    return RVIData{standardDeviation(highs), standardDeviation(lows), standardDeviation(closes)};
}

double calculateRVI(const RVIData& rviData) {
    double denominator = rviData.stdDevHigh + rviData.stdDevLow + rviData.stdDevClose;
    return (denominator == 0) ? 0 : rviData.stdDevClose / denominator;
}

//RSI
double calculateAverageGain(const std::deque<double>& gains, int period) {
    return std::accumulate(gains.begin(), gains.end(), 0.0) / period;
}

double calculateAverageLoss(const std::deque<double>& losses, int period) {
    return std::accumulate(losses.begin(), losses.end(), 0.0) / period;
}

double calculateRSI(const std::deque<double>& prices, int period) {
    std::deque<double> gains, losses;
    for (size_t i = 1; i < prices.size(); ++i) {
        double change = prices[i] - prices[i - 1];
        if (change > 0) gains.push_back(change);
        else losses.push_back(std::abs(change));
    }

    double avgGain = calculateAverageGain(gains, period);
    double avgLoss = calculateAverageLoss(losses, period);

    if (avgLoss == 0) return 100;
    double rs = avgGain / avgLoss;
    return 100 - (100 / (1 + rs));
}

double smiMax = 100;
double smiSignalLineMax = 98.0199;
double rviMax = .510228;
double rviSignalLineMax = .510228;
double rsiMax = 100;
double rsiSignalLineMax = 100;
double tempSpreadMax = .21;
double tempBookPressureMax = 0;
double tempTradingBalanceMax = .666619;
double adMax = 6.98231e+09;
double hvMax = 0.0245713;
double macdLineMax = -0.0111111;
double signalLineMax = -0.00924698;
double macdHistogramMax = -0.00888889;
double periodVolume2Max = 5.0353e+09;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string folderPath = "tickData/";
    std::string filePath = folderPath + argv[1];

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << filePath << std::endl;
        return 1;
    }

    std::string line;
    std::deque<TradeData> bidDeque, askDeque, highLowDeque, secondWindowDeque;
    std::deque<PriceData> priceDeque;
    std::deque<double> smiLineDeque; 
    std::deque<double> macdLineDeque; 
    std::deque<ResampledPriceData> resampledPriceDataDeque;
    std::deque<double> rviDeque; 
    std::deque<double> closingPrices; // for RSI
    std::deque<double> rsiDeque;

    std::getline(file, line); // Read the first line to skip the header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;

        std::getline(ss, cell, ',');
        std::tm dateTime = parseDateTime(cell);
        std::time_t currentTime = std::mktime(&dateTime);

        std::getline(ss, cell, ',');
        double close = std::stod(cell);
        priceDeque.push_back({currentTime, close});

        std::getline(ss, cell, ',');
        double tradeVolume = std::stod(cell);

        std::getline(ss, cell, ',');
        double totVol = std::stod(cell);
        highLowDeque.push_back({currentTime, close, totVol});
        secondWindowDeque.push_back({currentTime, 0, totVol});

        periodVolume += totVol;
        updateHighLowDeque(highLowDeque, currentTime, windowHigh, windowLow);
        updateSecondWindow(secondWindowDeque, currentTime);

        periodVolume = std::accumulate(highLowDeque.begin(), highLowDeque.end(), 0.0, [](double sum, const TradeData& trade) { return sum + trade.totVol; });
        periodVolume2 = std::accumulate(secondWindowDeque.begin(), secondWindowDeque.end(), 0.0, [](double sum, const TradeData& trade) { return sum + trade.totVol; });

        std::getline(ss, cell, ',');
        double bid = std::stod(cell);
        std::getline(ss, cell, ',');
        double ask = std::stod(cell);

        std::getline(ss, cell, ',');

        std::getline(ss, cell, ',');
        int bidSize = std::stoi(cell);

        std::getline(ss, cell, ',');
        int askSize = std::stoi(cell);

        for (int i = 0; i < 5; ++i) std::getline(ss, cell, ',');

        if (cell == "b") {
            bidDeque.push_back({currentTime, close});
            totalBidVolume += tradeVolume;
            updateDeque(bidDeque, currentTime, totalBidVolume);
        } else if (cell == "a") {
            askDeque.push_back({currentTime, close});
            totalAskVolume += tradeVolume;
            updateDeque(askDeque, currentTime, totalAskVolume);
        }

        double ad = calculateAccumulationDistribution(close, windowLow, windowHigh, periodVolume);
        double hv = calculateHistoricalVolatility(priceDeque);
        double shortTermEMA = calculateEMA(priceDeque, SHORT_TERM_DURATION);
        double longTermEMA = calculateEMA(priceDeque, LONG_TERM_DURATION);
        double macdLine = shortTermEMA - longTermEMA;
        macdLineDeque.push_back(macdLine);
        double signalLine = calculateEMAForMACD(macdLineDeque, SIGNAL_LINE_DURATION);
        double macdHistogram = macdLine - signalLine;

        std::map<std::time_t, ResampledPriceData> resampledPrices = resampleData(priceDeque);

        double smi = 0;
        double smiSignalLine = 0;
        double rvi = 0;
        double rviSignalLine = 0;
        double rsi = 0;
        double rsiSignalLine = 0;

        if (resampledPrices.size() >= SMI_LOOKBACK_PERIOD) {
            auto resampledPricesEnd = resampledPrices.end();
            double highestHigh = std::numeric_limits<double>::lowest();
            double lowestLow = std::numeric_limits<double>::max();

            std::advance(resampledPricesEnd, -SMI_LOOKBACK_PERIOD);
            for (auto it = resampledPricesEnd; it != resampledPrices.end(); ++it) {
                if (it->second.high > highestHigh) highestHigh = it->second.high;
                if (it->second.low < lowestLow) lowestLow = it->second.low;
            }

            double close = resampledPrices.rbegin()->second.close;
            double closeDiff = close - ((highestHigh + lowestLow) / 2);
            double rangeDiff = highestHigh - lowestLow;

            std::deque<double> closeDiffDeque;
            std::deque<double> rangeDiffDeque;
            closeDiffDeque.push_back(closeDiff);
            rangeDiffDeque.push_back(rangeDiff);

            double smoothedCloseDiff = calculateEMA(closeDiffDeque, SMI_SMOOTHING_PERIOD);
            double smoothedRangeDiff = calculateEMA(rangeDiffDeque, SMI_SMOOTHING_PERIOD);

            smi = 0;
            if (smoothedRangeDiff != 0) {
                smi = (smoothedCloseDiff / (smoothedRangeDiff / 2)) * 100;
            }
            smiLineDeque.push_back(smi);
            smiSignalLine = calculateEMA(smiLineDeque, SMI_SIGNAL_PERIOD);

            std::cout << "SMI: " << smi << ", "; 
            std::cout << "SMI Signal Line: " << smiSignalLine << ", ";

            if (smiMax < smi) smiMax = smi;
            if (smiSignalLineMax < smiSignalLine) smiSignalLineMax = smiSignalLine;
        }


        // RVI Calculation
        if (resampledPrices.size() >= RVI_PERIOD) {
            auto it = resampledPrices.end();
            std::deque<ResampledPriceData> recentData;
            
            std::advance(it, -RVI_PERIOD);

            for (; it != resampledPrices.end(); ++it) {
                recentData.push_back(it->second);
            }

            RVIData rviData = calculateRVIData(recentData);
            rvi = calculateRVI(rviData);
            rviDeque.push_back(rvi);

            if (rviDeque.size() > RVI_SIGNAL_PERIOD) {
                rviDeque.pop_front();
            }

            rviSignalLine = calculateEMA(rviDeque, RVI_SIGNAL_PERIOD);

            std::cout << "RVI: " << rvi << ", RVI Signal Line: " << rviSignalLine << ", ";

            if (rviMax < rvi) rviMax = rvi;
            if (rviSignalLineMax < rviSignalLine) rviSignalLineMax = rviSignalLine;

        }

        // RSI Calculation
        closingPrices.push_back(close);
        if (closingPrices.size() > RSI_PERIOD) {
            closingPrices.pop_front(); 
        }

        if (closingPrices.size() == RSI_PERIOD) {
            rsi = calculateRSI(closingPrices, RSI_PERIOD);
            rsiDeque.push_back(rsi);
            if (rsiDeque.size() > RSI_SIGNAL_PERIOD) {
                rsiDeque.pop_front();
            }
            rsiSignalLine = calculateEMA(rsiDeque, RSI_SIGNAL_PERIOD);

            std::cout << "RSI: " << rsi << ", RSI Signal Line: " << rsiSignalLine << ", ";

            if (rsiMax < rsi) rsiMax = rsi;
            if (rsiSignalLineMax < rsiSignalLine) rsiSignalLineMax = rsiSignalLine;

        }

        double tempSpread = calculateSpread(ask, bid);
        double tempBookPressure = calculateBookPressure(bidSize, askSize);
        double tempTradingBalance = calculateTradingBalance();
        
        std::cout << "Spread: " << tempSpread << ", ";
        std::cout << "Book Pressure: " << tempBookPressure << ", ";
        std::cout << "Trading Balance: " << tempTradingBalance << ", ";
        std::cout << "Accumulation/Distribution: " << ad << ", ";
        std::cout << "Historical Volatility: " << hv << ", ";
        std::cout << "MACD Line: " << macdLine << ", ";
        std::cout << "Signal Line: " << signalLine << ", ";
        std::cout << "MACD Histogram: " << macdHistogram << ", ";
        std::cout << "Period Volume: " << periodVolume2 << ", ";

        if (tempSpreadMax < tempSpread) tempSpreadMax = tempSpread;
        if (tempBookPressureMax > tempBookPressure) tempBookPressureMax = tempBookPressure;
        if (tempTradingBalanceMax < tempTradingBalance) tempTradingBalanceMax = tempTradingBalance;
        if (adMax < ad) adMax = ad;
        if (hvMax < hv) hvMax = hv;
        if (macdLineMax > macdLine) macdLineMax = macdLine;
        if (signalLineMax > signalLine) signalLineMax = signalLine;
        if (macdHistogramMax > macdHistogram) macdHistogramMax = macdHistogram;
        if (periodVolume2Max < periodVolume2) periodVolume2Max = periodVolume2;

        double wt = .07142857;

        
        double finalVal = smi/smiMax*wt + smiSignalLine/smiSignalLineMax*wt + rvi/rviMax*wt + rviSignalLine/rviSignalLineMax*wt + rsi/rsiMax*wt + rsiSignalLine/rsiSignalLineMax*wt + tempSpread/tempSpreadMax*wt + tempBookPressure/1*wt + tempTradingBalance/tempTradingBalanceMax*wt + ad/adMax*wt + hv/hvMax*wt + macdLine/macdLineMax*wt + signalLine/signalLineMax*wt + macdHistogram/macdHistogramMax*wt + periodVolume2/periodVolume2Max*wt;

        std::cout << "finalVal: " << finalVal << std::endl;
    }

    file.close();
    return 0;
}
