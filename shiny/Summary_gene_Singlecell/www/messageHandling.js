console.log('Shiny App loaded.');

// 函数检查Shiny是否已加载并且setInputValue是否可用
function isShinyReady() {
    return typeof Shiny !== 'undefined' && typeof Shiny.setInputValue === 'function';
}

var maxRetries = 10;
var retryCount = 0;
var initialDelay = 500;

function handleMessage(event) {
    if (isShinyReady()) {
        console.log('Shiny is ready and message received:', event.data);
        Shiny.setInputValue('data_from_html', event.data);
    } else if (retryCount < maxRetries) {
        var delay = initialDelay * (retryCount + 1);
        console.log(`Shiny is NOT ready yet. Retrying in ${delay}ms...`);
        setTimeout(function() {
            retryCount++;
            handleMessage(event);
        }, delay);
    } else {
        console.log('Max retries reached. Please check your Shiny app.');
    }
}

// 使用jQuery确保DOM已准备好
$(document).ready(function() {
    window.addEventListener('message', handleMessage);
});
