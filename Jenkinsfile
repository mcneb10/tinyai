pipeline {
    agent { label 'linux' }
    stages {
        stage('Setup') {
            steps {
                echo 'Updating...'
                sh 'apk update'
                echo '- Installing build tools'
                sh 'apk add --update alpine-sdk'
            }
        }
        stage('Build') {
            steps {
                echo 'Building...'
                sh 'make -C experiments/xor/'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing...'
                sh './experiments/xor/xor'
                sh './experiments/xor/xor_independent'
            }
        }
    }
}