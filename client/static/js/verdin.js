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
    console.log(results)
    resultTab.innerHTML = '<p class="text-danger">Return Json: ' + JSON.stringify(results) + '</p>'
}

function loadExample() {
    document.getElementById('variants').value = "1,25889628,1,25889638,DEL_3to5\n1,109792918,1,109792928,DEL_3to5\n"
}

function displayError(message) {
    resultTab.innerHTML = '<p class="text-danger">Error: ' + message + '</p>'
}
