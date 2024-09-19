console.log('Shiny App loaded.');

// 函数检查Shiny是否已加载并且setInputValue是否可用
function isShinyReady() {
    return typeof Shiny !== 'undefined' && typeof Shiny.setInputValue === 'function';
}

function handleMessage(event) {
    setTimeout(function() {
        if (isShinyReady()) {
            console.log('Shiny is ready and message received:', event.data);
            Shiny.setInputValue('data_from_html', event.data);
        } else {
            console.log('Shiny is NOT ready yet.');
        }
    }, 500);  // 等待 500 毫秒
}

// 使用jQuery确保DOM已准备好
$(document).ready(function() {
    window.addEventListener('message', handleMessage);
});

