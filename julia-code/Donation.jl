import HTTP
import JSON

using Dates

urlg = "https://ssl.gongyi.qq.com/cgi-bin/gywcom_proj_trans_query_ch?pid=218649&type=detail"
urlf = "https://ssl.gongyi.qq.com/cgi-bin/gywcom_proj_trans_query_hb?pid=218649&lst="
outpath = "C:\\Users\\zhangyutong926\\Desktop\\res\\res.csv"
f = open(outpath, "w")

lst = ""
while true
    url = lst == "" ? urlg : urlf * lst
    println(url)
    req = HTTP.request("GET", url)
    json = String(req.body)
    dict = JSON.parse(json)
    code = get(dict, "code", -1)
    msg = get(dict, "msg", "")
    optime = get(dict, "op_time", -1)
    caer = get(dict["data"], "caer", -1)
    global lst = dict["data"]["lst"]
    if get(dict["data"], "trans", Dict()) == Dict()
        break
    else
        for i in dict["data"]["trans"]
            money = i["m"] / 100.
            timeElapsed = i["t"]
            realtime = optime - timeElapsed
            name = i["n"]
            header = i["h"]
            item = "$code, \"$msg\", $optime, $money, $realtime, \"$name\", \"$header\", $caer"
            write(f, item * "\n")
        end
    end
    sleep(.1)
end

close(f)
