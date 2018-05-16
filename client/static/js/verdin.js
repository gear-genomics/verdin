/* global c3 */

$('#mainTab a').on('click', function (e) {
    e.preventDefault()
    $(this).tab('show')
})

var resultLink = document.getElementById('link-results')
var resultTab = document.getElementById('result-tab')

var submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function () {
    resultLink.click()
    run()
})

var exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', function () {
    loadExample()
})

var varsInput = document.getElementById('variants')

var spinnerHtml = '<i class="fas fa-spinner fa-2x spinner"></i>'

function run() {
    resultTab.innerHTML = spinnerHtml

    // Parse variants as json
    var variants = varsInput.value.split('\n').filter(function (line) { return line !== ""; })
    var arr = []
    for (var i =0; i < variants.length; ++i) {
	fields = variants[i].split(',')
	var jsdict = { "chr1": fields[0], "pos1": parseInt(fields[1]), "chr2": fields[2], "pos2": parseInt(fields[3]), "type": fields[4] }
	arr.push(jsdict)
    }

    // Send json post request
    var req = new XMLHttpRequest()
    var url = "http://0.0.0.0:3300/primers"
    var data = JSON.stringify(arr)
    req.addEventListener('load', displayResults)
    req.open('POST', url, true)
    req.setRequestHeader("Content-type", "application/json")
    req.send(data)

    // Get request
    //var req = new XMLHttpRequest()
    //req.addEventListener('load', displayResults)
    //req.open('GET', 'http://0.0.0.0:3300/primers?chr1=1&pos1=7878&chr2=7&pos2=56889&type=BND_3to5')
    //req.send()
}

function displayResults() {
    var results = JSON.parse(this.responseText)
    //console.log(results)
    //resultTab.innerHTML = '<p class="text-danger">Return Json: ' + JSON.stringify(results) + '</p>'
    var str = '<table class="table">'
    str += '<thead><tr>'
    str += '<th scope="col">chr1</th>'
    str += '<th scope="col">pos1</th>'
    str += '<th scope="col">chr2</th>'
    str += '<th scope="col">pos2</th>'
    str += '<th scope="col">type</th>'
    str += '<th scope="col">Primer1Seq</th>'
    str += '<th scope="col">Primer1Chrom</th>'
    str += '<th scope="col">Primer1Pos</th>'
    str += '<th scope="col">Primer1Ori</th>'
    str += '<th scope="col">Primer1Tm</th>'
    str += '<th scope="col">Primer2Seq</th>'
    str += '<th scope="col">Primer2Chrom</th>'
    str += '<th scope="col">Primer2Pos</th>'
    str += '<th scope="col">Primer2Ori</th>'
    str += '<th scope="col">Primer2Tm</th>'
    str += '</tr></thead>'
    str += '<tbody>'
    for (var i = 0; i < results.length; ++i) {
	str += '<tr>'
	str += '<td>' + results[i]['chr1'] + '</td>'
	str += '<td>' + results[i]['pos1'] + '</td>'
	str += '<td>' + results[i]['chr2'] + '</td>'
	str += '<td>' + results[i]['pos2'] + '</td>'
	str += '<td>' + results[i]['type'] + '</td>'
	str += '<td>' + results[i]['Primer1Seq'] + '</td>'
	str += '<td>' + results[i]['Primer1Chrom'] + '</td>'
	str += '<td>' + results[i]['Primer1Pos'] + '</td>'
	str += '<td>' + results[i]['Primer1Ori'] + '</td>'
	str += '<td>' + results[i]['Primer1Tm'] + '</td>'
	str += '<td>' + results[i]['Primer2Seq'] + '</td>'
	str += '<td>' + results[i]['Primer2Chrom'] + '</td>'
	str += '<td>' + results[i]['Primer2Pos'] + '</td>'
	str += '<td>' + results[i]['Primer2Ori'] + '</td>'
	str += '<td>' + results[i]['Primer2Tm'] + '</td>'	
	str += '</tr>'
    }
    str += '</tbody>'
    str += '</table>'
    resultTab.innerHTML = str
}

function loadExample() {
    document.getElementById('variants').value = "1,201295850,1,202301202,DEL\n2,134118077,2,137059597,INV_5to5\n2,205249217,20,57401265,BND_3to5\n10,115037918,10,115047289,DUP\n11,113564095,11,113564096,SNV\n11,114466285,11,116054519,INV_3to3\n12,71830272,4,8664772,BND_5to3\n13,35642391,6,46761397,BND_5to5\n13,35651962,6,46768987,BND_3to3\n17,65888129,17,65888129,INS\n21,38797660,21,38801104,DEL\n"
}

function displayError(message) {
    resultTab.innerHTML = '<p class="text-danger">Error: ' + message + '</p>'
}
