pipeline {
    agent { 
        docker {
            image 'ubuntu'
            args '-u root:sudo'
        }
    }
    stages {
        stage('Setup') {
            def dockerHome = tool 'MyDocker'
            env.PATH = "${dockerHome}/bin:${env.PATH}"
            steps {
                echo 'Updating...'
                sh 'apt update'
                echo '- Installing build tools'
                sh 'apt install build-essential'
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