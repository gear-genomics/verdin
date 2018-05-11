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
    req.open('GET', '/chr1/7878/chr7/56889/BND_3to5')
    req.send()
}

function displayResults() {
    var results = JSON.parse(this.responseText)
    resultTab.innerHTML = '<p class="text-danger">Return Json: ' + JSON.stringify(results) + '</p>'
}

function loadExample() {
    resultTab.innerHTML = spinnerHtml
    var req = new XMLHttpRequest()
    req.addEventListener('load', displayResults)
    req.open('GET', '/chr1/7878/chr7/56889/BND_3to5')
    req.send()
}

function displayError(message) {
    resultTab.innerHTML = '<p class="text-danger">Error: ' + message + '</p>'
}
