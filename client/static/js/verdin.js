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
    resultLink.click()
    loadExample()
})

var varsInput = document.getElementById('variants')

var spinnerHtml = '<i class="fas fa-spinner fa-2x spinner"></i>'

function run() {
    resultTab.innerHTML = spinnerHtml
    //ToDo: Parse input text fields and generate URL
    var req = new XMLHttpRequest()
    req.addEventListener('load', displayResults)
    req.open('GET', 'http://localhost:3300/primers?build=hg19&chr1=chr1&pos1=7878&chr2=chr7&pos2=56889&svtype=BND_3to5')
    req.send()
}

function displayResults() {
    var results = JSON.parse(this.responseText)
    console.log(results)
    resultTab.innerHTML = '<p class="text-danger">Return Json: ' + JSON.stringify(results) + '</p>'
}


var exampleVariants = [
    {"chr1": "chr7", "pos1": 20300000, "chr2": "chr7", "pos2": 20600000, "type": "DEL_3to5"},
    {"chr1": "chr7", "pos1": 40300000, "chr2": "chr7", "pos2": 40600000, "type": "DEL_3to5"}
]


function loadExample() {
    resultTab.innerHTML = spinnerHtml
    var req = new XMLHttpRequest()
    var url = "http://0.0.0.0:3300/primers"
    var data = JSON.stringify(exampleVariants)
    req.addEventListener('load', displayResults)
    req.open('POST', url, true)
    req.setRequestHeader("Content-type", "application/json")
    req.send(data)
}

function displayError(message) {
    resultTab.innerHTML = '<p class="text-danger">Error: ' + message + '</p>'
}
